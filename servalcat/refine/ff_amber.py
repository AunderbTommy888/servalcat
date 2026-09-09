"""
AMBER force-field prior for Servalcat refinement (Phase 1).

This module is intentionally minimal:
- Optional dependency on OpenMM
- Energy + Cartesian gradient contribution only
- XYZ refinement only (no ADP/occ coupling)
- Alternative conformations: only the first conformer of each residue is
  passed to OpenMM; the remaining altloc atoms receive zero force-field gradient.
- Histidine protonation state can be chosen (HIP/HIE/HID) by excluding HD1/HE2
  from the force-field model.
- Hessian is approximated by its diagonal, estimated per atom from the bonded
  (bond + angle) terms in the Gauss-Newton sense, with a constant floor.
"""
from __future__ import absolute_import, division, print_function, generators
import os
import tempfile
import numpy
import scipy.sparse
import gemmi
from servalcat import ext
from servalcat.utils import logger
from servalcat.refine import ff_ligand

HIS_STATES = ("HIP", "HIE", "HID")
HESSIAN_MODES = ("bonded", "const")
# The angle gradient carries a 1/sin(theta) factor. It stays finite only away from a linear
# arrangement, so the factor is floored: a (near-)linear angle would otherwise contribute a block
# with an arbitrarily large norm, which destroys the conditioning of the whole normal matrix and,
# in floating point, can break the positive semi-definiteness that holds algebraically.
MIN_SIN_ANGLE = 0.08715574  # sin(5 degrees)


def xyz_grad_to_param_grad(grad_xyz, atom_to_param, n_params):
    """Map per-atom xyz gradients to refinement parameter vector.

    Parameters
    ----------
    grad_xyz : (n_atom, 3) ndarray
        dE/dx in energy per angstrom.
    atom_to_param : sequence[int]
        Mapping from atom index to xyz-block index in refinement params.
    n_params : int
        Total number of refinement parameters.
    """
    g = numpy.zeros(n_params, dtype=float)
    if len(atom_to_param) != grad_xyz.shape[0]:
        raise RuntimeError("atom_to_param and grad_xyz sizes are inconsistent")
    for i_atom, pidx in enumerate(atom_to_param):
        if pidx < 0:
            continue
        j = int(pidx) * 3
        if j + 2 >= n_params:
            raise RuntimeError("xyz parameter mapping is out of range")
        g[j:j+3] = grad_xyz[i_atom]
    return g


def select_first_conformer(model):
    """Select atoms that form the first alternative conformation.

    Returns a list of bool flags in the order of ``model.all()``.
    Rules (same spirit as gemmi.Structure.remove_alternative_conformations):
    - atoms without altloc are always kept
    - among altloc atoms of a residue, the first atom of each name is kept
    - a residue sharing seqid with the preceding residue in the chain
      (microheterogeneity) is treated as an alternative residue and dropped
    """
    keep = []
    for chain in model:
        prev_seqid = None
        for res in chain:
            alt_residue = (prev_seqid is not None and res.seqid == prev_seqid)
            prev_seqid = res.seqid
            seen = set()
            for atom in res:
                if alt_residue:
                    keep.append(False)
                elif atom.altloc == '\0':
                    keep.append(True)
                elif atom.name in seen:
                    keep.append(False)
                else:
                    seen.add(atom.name)
                    keep.append(True)
    return keep


def select_his_protonation(model, his_state):
    """Flags (order of ``model.all()``) of atoms to exclude from the force field
    so that histidines match the requested AMBER protonation state.

    HIP: keep HD1 and HE2 (as generated from the monomer library)
    HIE: drop HD1 (proton on NE2 only)
    HID: drop HE2 (proton on ND1 only)
    """
    if his_state not in HIS_STATES:
        raise RuntimeError("Unknown histidine state: {} (choose from {})".format(his_state, HIS_STATES))
    drop_name = {"HIP": None, "HIE": "HD1", "HID": "HE2"}[his_state]
    exclude = []
    for chain in model:
        for res in chain:
            for atom in res:
                exclude.append(drop_name is not None and res.name == "HIS" and atom.name == drop_name)
    return exclude


def _keep_positive(term, k_index):
    """Filter a (i, j, [k], force_constant) tuple of arrays down to the terms with a positive
    force constant. Returns (filtered, n_dropped)."""
    k = term[k_index]
    keep = k > 0
    n_dropped = int(numpy.count_nonzero(~keep))
    if not n_dropped:
        return term, 0
    return tuple(a[keep] for a in term), n_dropped


def angle_gradients(pos, ai, aj, ak):
    """Gradients of the angle i-j-k (j central) with respect to the three atoms.

    grad_i theta = (cos t * a_hat - b_hat) / (|a| sin t),  a = r_i - r_j
    grad_k theta = (cos t * b_hat - a_hat) / (|b| sin t),  b = r_k - r_j
    grad_j theta = -(grad_i theta + grad_k theta)

    Returns (grad_i, grad_j, grad_k, n_near_linear); sin t is floored at MIN_SIN_ANGLE.
    """
    a = pos[ai] - pos[aj]
    b = pos[ak] - pos[aj]
    la = numpy.linalg.norm(a, axis=1)
    lb = numpy.linalg.norm(b, axis=1)
    ah = a / la[:, None]
    bh = b / lb[:, None]
    cos_t = numpy.clip(numpy.sum(ah * bh, axis=1), -1., 1.)
    sin_t = numpy.sqrt(numpy.maximum(1. - cos_t**2, 0.))
    n_near_linear = int(numpy.count_nonzero(sin_t < MIN_SIN_ANGLE))
    sin_t = numpy.maximum(sin_t, MIN_SIN_ANGLE)
    gi = (cos_t[:, None] * ah - bh) / (la * sin_t)[:, None]
    gk = (cos_t[:, None] * bh - ah) / (lb * sin_t)[:, None]
    return gi, -(gi + gk), gk, n_near_linear


def bonded_hessian_blocks(pos, bonds, angles, atom_to_param):
    """Gauss-Newton blocks of the bonded Hessian as (param_idx_a, param_idx_b, (M,3,3)) tuples.

    Each term is E = k/2 q(x)^2, so its Hessian is k (grad q)(grad q)^T + k q grad^2 q; the second
    part is dropped (Gauss-Newton). What is left is a sum of k * outer(g, g) with k >= 0, which is
    positive semi-definite by construction - the property the normal matrix relies on.

    Bond  (i,j): H_ii = H_jj = k u u^T, H_ij = -k u u^T with u the bond direction
    Angle (i,j,k): H_mn = k_theta (grad_m theta)(grad_n theta)^T for the six upper-triangle pairs
    """
    a2p = numpy.asarray(atom_to_param, dtype=int)
    out = []
    bi, bj, bk = bonds
    if len(bi):
        v = pos[bi] - pos[bj]
        u = v / numpy.linalg.norm(v, axis=1)[:, None]
        H = bk[:, None, None] * u[:, :, None] * u[:, None, :]
        out += [(a2p[bi], a2p[bi], H), (a2p[bj], a2p[bj], H), (a2p[bi], a2p[bj], -H)]
    ai, aj, ak, akth = angles
    if len(ai):
        gi, gj, gk, _ = angle_gradients(pos, ai, aj, ak)
        ps, gs = (a2p[ai], a2p[aj], a2p[ak]), (gi, gj, gk)
        for m in range(3):
            for l in range(m, 3):
                out.append((ps[m], ps[l], akth[:, None, None] * gs[m][:, :, None] * gs[l][:, None, :]))
    return out


def assemble_bonded_hessian(pos, bonds, angles, atom_to_param, n_params, diag_floor=0.):
    """Sparse Gauss-Newton Hessian of the bonded terms on the parameter vector.

    Positive semi-definite by construction (see bonded_hessian_blocks); diag_floor adds a
    non-negative diagonal, which cannot make it indefinite and lifts the rigid-body null modes.
    """
    rows, cols, vals = [], [], []
    for pa, pb, H in bonded_hessian_blocks(pos, bonds, angles, atom_to_param):
        same = pa is pb or numpy.array_equal(pa, pb)
        ok = (pa >= 0) & (pb >= 0)
        if not ok.any():
            continue
        pa_, pb_, H_ = pa[ok], pb[ok], H[ok]
        k = numpy.arange(3)
        r = (pa_[:, None, None] * 3 + k[None, :, None]) * numpy.ones((1, 1, 3), dtype=int)
        c = (pb_[:, None, None] * 3 + k[None, None, :]) * numpy.ones((1, 3, 1), dtype=int)
        rows.append(r.ravel()); cols.append(c.ravel()); vals.append(H_.ravel())
        if not same:
            # the transposed block: rows/cols are already swapped, so the values stay as they are
            # (A[3*pb+b, 3*pa+a] = H[a, b] makes the assembled matrix symmetric)
            rows.append(c.ravel()); cols.append(r.ravel()); vals.append(H_.ravel())
    if not rows:
        A = scipy.sparse.csr_matrix((n_params, n_params))
    else:
        A = scipy.sparse.coo_matrix((numpy.concatenate(vals),
                                    (numpy.concatenate(rows), numpy.concatenate(cols))),
                                   shape=(n_params, n_params)).tocsr()
    if diag_floor:
        touched = numpy.zeros(n_params, dtype=bool)
        a2p = numpy.asarray(atom_to_param, dtype=int)
        for p in a2p[a2p >= 0]:
            touched[p*3:p*3+3] = True
        add = numpy.where(touched, numpy.maximum(diag_floor - A.diagonal(), 0.), 0.)
        A = A + scipy.sparse.diags(add, format="csr")
    return A


def bonded_hessian_diag(pos, bonds, angles):
    """Gauss-Newton diagonal of the Hessian from harmonic bond and angle terms.

    pos    : (n, 3) positions in angstrom
    bonds  : (i, j, k) arrays; k in kJ/mol/A^2
    angles : (i, j, k, kth) arrays; j is the central atom; kth in kJ/mol/rad^2

    Bond  E = k/2 (r - r0)^2      -> d2E/dx_{i,a}^2 ~ k * u_a^2, u = (r_i - r_j)/|r_i - r_j|
    Angle E = kth/2 (t - t0)^2    -> d2E/dx_{m,a}^2 ~ kth * (dt/dx_{m,a})^2
      grad_i t = (cos t * a_hat - b_hat) / (|a| sin t),  a = r_i - r_j
      grad_k t = (cos t * b_hat - a_hat) / (|b| sin t),  b = r_k - r_j
      grad_j t = -(grad_i t + grad_k t)
    Returns (n, 3) array in kJ/mol/A^2.
    """
    d = numpy.zeros_like(pos)
    bi, bj, bk = bonds
    if len(bi):
        v = pos[bi] - pos[bj]
        r = numpy.linalg.norm(v, axis=1)
        u2 = (v / r[:, None])**2
        numpy.add.at(d, bi, bk[:, None] * u2)
        numpy.add.at(d, bj, bk[:, None] * u2)
    ai, aj, ak, akth = angles
    if len(ai):
        gi, gj, gk, _ = angle_gradients(pos, ai, aj, ak)
        numpy.add.at(d, ai, akth[:, None] * gi**2)
        numpy.add.at(d, ak, akth[:, None] * gk**2)
        numpy.add.at(d, aj, akth[:, None] * gj**2)
    return d


class AmberFFPrior:
    def __init__(self, st, refine_params,
                 forcefield_files=None,
                 nonbonded_method="NoCutoff",
                 platform_name="Reference",
                 hessian_diag=1000.,
                 hessian_mode="bonded",
                 his_state="HIP",
                 monlib=None,
                 ligand_ff=ff_ligand.DEFAULT_LIGAND_FF,
                 ligand_charge=ff_ligand.DEFAULT_CHARGE_METHOD,
                 ligand_smiles=None):
        self.st = st
        self.monlib = monlib
        self.ligand_ff = ligand_ff
        self.ligand_charge = ligand_charge
        self.ligand_smiles = dict(ligand_smiles or {})
        self.ligand_info = []
        self.n_disulfides = 0
        self.geom_weights_full = None  # per-atom geometry weights before replace_geom_restraints()
        self.refine_params = refine_params
        self.hessian_diag = float(hessian_diag)  # constant value (const) or floor (bonded), kJ/mol/A^2
        if self.hessian_diag < 0:
            raise RuntimeError("hessian_diag must not be negative (got {}): a negative diagonal would "
                               "make the normal matrix indefinite".format(hessian_diag))
        self._warned_near_linear = False
        if hessian_mode not in HESSIAN_MODES:
            raise RuntimeError("Unknown hessian mode: {} (choose from {})".format(hessian_mode, HESSIAN_MODES))
        self.hessian_mode = hessian_mode
        self.his_state = his_state
        if nonbonded_method not in ("NoCutoff", "CutoffNonPeriodic"):
            raise RuntimeError("Unsupported nonbonded method for AMBER prior: {} (periodic methods are not supported; "
                               "use NoCutoff or CutoffNonPeriodic)".format(nonbonded_method))
        self._openmm = None
        self._unit = None
        self._context = None
        self._integrator = None
        self._n_atoms = len(self.refine_params.atoms)
        self._ordered_atoms = [cra.atom for cra in self.st[0].all()]

        # Atoms passed to OpenMM: first conformer only, minus histidine protons to drop
        keep_alt = select_first_conformer(self.st[0])
        drop_his = select_his_protonation(self.st[0], his_state)
        if len(keep_alt) != len(self._ordered_atoms) or len(drop_his) != len(self._ordered_atoms):
            raise RuntimeError("Internal error: atom traversal order mismatch")
        self._keep = [k and not d for k, d in zip(keep_alt, drop_his)]
        self._ff_atoms = [a for a, k in zip(self._ordered_atoms, self._keep) if k]
        self.n_ff_atoms = len(self._ff_atoms)
        self.n_excluded_atoms = len(self._ordered_atoms) - self.n_ff_atoms
        n_alt = keep_alt.count(False)
        n_his = sum(1 for k, d in zip(keep_alt, drop_his) if k and d)
        if n_alt > 0:
            logger.writeln("AMBER prior: {} atom(s) in alternative conformations are excluded "
                           "from the force field (first conformer only)".format(n_alt))
        if n_his > 0:
            logger.writeln("AMBER prior: {} histidine proton(s) excluded from the force field "
                           "to use the {} state".format(n_his, his_state))
        self._check_c_termini()

        if forcefield_files is None:
            forcefield_files = ["amber14-all.xml", "amber14/tip3p.xml"]
        self.forcefield_files = list(forcefield_files)
        self.nonbonded_method = nonbonded_method
        self.platform_name = platform_name

        self._setup_openmm_system()

    def _check_c_termini(self):
        """Warn when the last polymer residue of a chain lacks OXT: OpenMM then treats it as an
        internal residue (with ignoreExternalBonds=True) instead of a C-terminal template."""
        for chain in self.st[0]:
            poly = [r for r in chain if r.het_flag != 'H' and r.name != "HOH"]
            if not poly:
                continue
            last = poly[-1]
            if last.find_atom("OXT", '*') is None and last.find_atom("C", '*') is not None:
                logger.writeln("AMBER prior: WARNING: {}/{} {} has no OXT; treated as an internal residue "
                               "with a dangling C (no C-terminal template/charges)".format(chain.name, last.name, last.seqid))

    def _make_ff_structure(self):
        """Clone the structure keeping only the atoms selected in self._keep."""
        st_ff = self.st.clone()
        i = 0
        for chain in st_ff[0]:
            for res in chain:
                n = len(res)
                flags = self._keep[i:i+n]
                i += n
                for j in reversed(range(n)):
                    if flags[j]:
                        res[j].altloc = '\0'
                    else:
                        del res[j]
        st_ff.remove_empty_chains()
        return st_ff

    def _setup_openmm_system(self):
        try:
            try:
                import openmm
                from openmm import app, unit
            except ImportError:
                from simtk import openmm
                from simtk.openmm import app
                from simtk import unit
        except ImportError as e:
            raise RuntimeError("OpenMM is required for --amber_enable") from e

        self._openmm = openmm
        self._unit = unit

        st_ff = self._make_ff_structure()
        with tempfile.TemporaryDirectory(prefix="servalcat_amber_") as dtmp:
            pdb_path = os.path.join(dtmp, "amber_input.pdb")
            st_ff.write_pdb(pdb_path)
            pdb = app.PDBFile(pdb_path)
            # The unit cell of an SPA model is the map box; never let it act as a periodic cell.
            pdb.topology.setPeriodicBoxVectors(None)
            # Disulfide-bonded cysteines (SG-SG bonds are created by PDBFile from distances): with
            # ignoreExternalBonds=True a CYS without HG matches both CYM and CYX, so pin them to CYX.
            residue_templates = {}
            for a1, a2 in pdb.topology.bonds():
                if (a1.residue is not a2.residue and a1.name == "SG" and a2.name == "SG"
                        and a1.residue.name == "CYS" and a2.residue.name == "CYS"):
                    residue_templates[a1.residue] = "CYX"
                    residue_templates[a2.residue] = "CYX"
            self.n_disulfides = len(residue_templates) // 2
            if self.n_disulfides:
                logger.writeln("AMBER prior: {} disulfide bond(s) -> CYX templates".format(self.n_disulfides))

            ff = app.ForceField(*self.forcefield_files)
            if not hasattr(app, self.nonbonded_method):
                raise RuntimeError("Unknown OpenMM nonbonded method: {}".format(self.nonbonded_method))
            nb_method = getattr(app, self.nonbonded_method)
            if self.ligand_ff and str(self.ligand_ff).lower() != "none":
                if ff_ligand.find_ligand_residues(ff, pdb.topology):
                    if self.monlib is None:
                        raise RuntimeError("residues without AMBER templates found but no monomer library was given "
                                           "to AmberFFPrior (needed for ligand force-field assignment)")
                    self.ligand_info = ff_ligand.register_ligand_templates(
                        ff, pdb.topology, self.monlib, ligand_ff=self.ligand_ff,
                        charge_method=self.ligand_charge, smiles_overrides=self.ligand_smiles)
                    for name, n, q, lff in self.ligand_info:
                        logger.writeln("AMBER prior: {} parameterised with {} ({} atoms, formal charge {:+d}, {} charges)".format(
                            name, lff, n, q, self.ligand_charge))
            try:
                system = ff.createSystem(pdb.topology,
                                         nonbondedMethod=nb_method,
                                         constraints=None,
                                         rigidWater=False,
                                         residueTemplates=residue_templates,
                                         ignoreExternalBonds=True)
            except Exception as e:
                raise RuntimeError(
                    "OpenMM failed to build an AMBER system. "
                    "This usually means the model contains residues without force-field templates "
                    "(e.g. ligands); enable ligand parameterisation with --amber_ligand_ff (needs rdkit, "
                    "openff-toolkit and openmmforcefields) or supply parameter XML files via --amber_forcefield.\n"
                    "Original error: {}".format(e)
                ) from e

            self._integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
            if self.platform_name:
                platform = openmm.Platform.getPlatformByName(self.platform_name)
                self._context = openmm.Context(system, self._integrator, platform)
            else:
                self._context = openmm.Context(system, self._integrator)

            self._topology_atoms = list(pdb.topology.atoms())
            self._topology_n_atoms = len(self._topology_atoms)
            self._extract_bonded_terms(system)

        if self._topology_n_atoms != self.n_ff_atoms:
            raise RuntimeError(
                "OpenMM topology atom count ({}) does not match structure atom count ({})".format(
                    self._topology_n_atoms, self.n_ff_atoms
                )
            )

        # RefineParams.atom_to_param is indexed by atom.serial - 1
        atom_to_param_ref = list(self.refine_params.atom_to_param(ext.RefineParams.Type.X))
        self._atom_to_param_x = []
        for atom in self._ff_atoms:
            idx = atom.serial - 1
            if idx < 0 or idx >= len(atom_to_param_ref):
                raise RuntimeError("Atom serial {} is out of range for refinement parameters; "
                                   "atom serials must be renumbered 1..N before AmberFFPrior".format(atom.serial))
            self._atom_to_param_x.append(atom_to_param_ref[idx])

    def _extract_bonded_terms(self, system):
        """Collect harmonic bond/angle parameters for the diagonal Hessian estimate."""
        openmm, unit = self._openmm, self._unit
        bi, bj, bk = [], [], []
        ai, aj, ak, akth = [], [], [], []
        for force in system.getForces():
            if isinstance(force, openmm.HarmonicBondForce):
                for n in range(force.getNumBonds()):
                    i, j, r0, k = force.getBondParameters(n)
                    bi.append(i); bj.append(j)
                    # kJ/mol/nm^2 -> kJ/mol/A^2
                    bk.append(k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer**2) * 0.01)
            elif isinstance(force, openmm.HarmonicAngleForce):
                for n in range(force.getNumAngles()):
                    i, j, k, t0, kth = force.getAngleParameters(n)
                    ai.append(i); aj.append(j); ak.append(k)
                    akth.append(kth.value_in_unit(unit.kilojoule_per_mole / unit.radian**2))
        bonds = (numpy.array(bi, dtype=int), numpy.array(bj, dtype=int), numpy.array(bk, dtype=float))
        angles = (numpy.array(ai, dtype=int), numpy.array(aj, dtype=int), numpy.array(ak, dtype=int),
                  numpy.array(akth, dtype=float))
        # A Gauss-Newton block k*outer(g, g) is positive semi-definite only for k >= 0. Force fields
        # should never give a negative force constant, but drop such terms rather than let them make
        # the normal matrix indefinite. Zero constants contribute nothing and are dropped too.
        self._bonds, n_bad_b = _keep_positive(bonds, 2)
        self._angles, n_bad_a = _keep_positive(angles, 3)
        if n_bad_b or n_bad_a:
            logger.writeln("AMBER prior: WARNING: dropped {} bond and {} angle term(s) with a "
                           "non-positive force constant from the Hessian estimate".format(n_bad_b, n_bad_a))

    def _positions_angstrom(self):
        return numpy.array([[a.pos.x, a.pos.y, a.pos.z] for a in self._ff_atoms], dtype=float)

    def _positions_nm(self):
        # OpenMM expects nm; gemmi positions are angstrom.
        out = []
        for atom in self._ff_atoms:
            p = atom.pos
            out.append(self._openmm.Vec3(p.x * 0.1, p.y * 0.1, p.z * 0.1))
        return out

    def hessian_diag_atoms(self):
        """Per-atom (n_ff_atoms, 3) diagonal Hessian estimate in kJ/mol/A^2 (unweighted)."""
        if self.hessian_mode == "const":
            return numpy.full((self.n_ff_atoms, 3), self.hessian_diag)
        d = bonded_hessian_diag(self._positions_angstrom(), self._bonds, self._angles)
        return numpy.maximum(d, self.hessian_diag)

    def hessian_diag_vector(self):
        """Diagonal Hessian estimate mapped onto the refinement parameter vector (unweighted).
        Parameters of atoms outside the force field (excluded altloc/protons, ADP, occ) get 0."""
        n_params = self.refine_params.n_params()
        d = numpy.zeros(n_params, dtype=float)
        da = self.hessian_diag_atoms()
        for i_atom, pidx in enumerate(self._atom_to_param_x):
            if pidx < 0:
                continue
            j = int(pidx) * 3
            d[j:j+3] = da[i_atom]
        return d

    def hessian_blocks(self):
        """Gauss-Newton blocks of the bonded Hessian at the current coordinates."""
        return bonded_hessian_blocks(self._positions_angstrom(), self._bonds, self._angles,
                                     self._atom_to_param_x)

    def hessian_matrix(self):
        """Sparse Hessian estimate on the refinement parameter vector, including the off-diagonal
        blocks that couple bonded atoms. The diagonal is identical to hessian_diag_vector().

        Positive semi-definite by construction: a sum of k*outer(g, g) with k >= 0 (non-positive
        force constants are dropped in _extract_bonded_terms) plus a non-negative diagonal floor.
        """
        if self.hessian_mode != "bonded":
            raise RuntimeError("off-diagonal Hessian is only available with hessian_mode='bonded' "
                               "(got '{}')".format(self.hessian_mode))
        pos = self._positions_angstrom()
        self._warn_near_linear_angles(pos)
        A = assemble_bonded_hessian(pos, self._bonds, self._angles, self._atom_to_param_x,
                                    self.refine_params.n_params())
        # same floor as hessian_diag_vector(): the bonded Hessian is rank-deficient (rigid-body
        # modes), so without it the total matrix can be singular. A non-negative diagonal cannot
        # make a positive semi-definite matrix indefinite.
        dv = self.hessian_diag_vector()
        return A + scipy.sparse.diags(numpy.maximum(dv - A.diagonal(), 0.), format="csr")

    def _warn_near_linear_angles(self, pos):
        """Report once if any angle is close enough to linear for the 1/sin(theta) floor to bite."""
        ai, aj, ak, _ = self._angles
        if not len(ai) or self._warned_near_linear:
            return
        _, _, _, n = angle_gradients(pos, ai, aj, ak)
        if n:
            self._warned_near_linear = True
            logger.writeln("AMBER prior: WARNING: {} angle(s) within {:.1f} deg of linear; the 1/sin(theta) "
                           "factor of the Hessian estimate is clamped there".format(
                               n, numpy.degrees(numpy.arcsin(MIN_SIN_ANGLE))))

    def replace_geom_restraints(self, ratio):
        """Hand the classical geometry restraints over to the force field.

        Servalcat weights every restraint by the mean of the per-atom weights of the atoms it
        involves (RefineParams::find_geom_weight), so multiplying the per-atom weight of the atoms
        this force field covers by (1 - ratio) scales those restraints down:
        ratio = 1 removes them entirely and AMBER alone provides the chemistry, ratio = 0 leaves
        them as they are and AMBER is an additional term.

        Atoms outside the force field keep their full classical restraints - alternative conformers
        past the first, and histidine protons dropped for HIE/HID. A restraint that spans both kinds
        of atom keeps the average of the two weights, which is what makes this graceful.

        Not covered by AMBER, and therefore lost in proportion to ratio (the per-atom weight applies
        to every restraint type):
        - chirality: AMBER has no explicit sp3 chirality term (its impropers only keep sp2 planarity)
        - NCS and stacking restraints
        - VDW repulsion against symmetry and NCS copies: OpenMM only ever sees the asymmetric unit
        ADP, occupancy, jelly-body and external restraints are unaffected (they do not go through
        the per-atom geometry weight).

        Returns the number of atoms whose weight was scaled.
        """
        if not 0. <= ratio <= 1.:
            raise RuntimeError("geometry replacement ratio must be within [0, 1] (got {})".format(ratio))
        if ratio == 0.:
            return 0
        w = self.refine_params.geom_weights  # per-atom view, indexed by atom.serial - 1
        # keep the original weights: zero-weight restraints drop out of the geometry report, so the
        # reporting pass restores them to keep bond/angle rmsZ available as a diagnostic
        self.geom_weights_full = numpy.array(w, copy=True)
        serials = numpy.array([a.serial for a in self._ff_atoms], dtype=int) - 1
        if serials.size and (serials.min() < 0 or serials.max() >= len(w)):
            raise RuntimeError("atom serials are out of range for the geometry weight vector; "
                               "they must be renumbered 1..N before replace_geom_restraints()")
        w[serials] *= (1. - ratio)  # multiply, so a local_geom_weights setting is not discarded
        logger.writeln("AMBER prior: classical geometry restraints scaled by {:.3g} for {} atom(s) "
                       "covered by the force field{}".format(
                           1. - ratio, serials.size,
                           "; {} atom(s) outside it keep full restraints".format(self.n_excluded_atoms)
                           if self.n_excluded_atoms else ""))
        if ratio > 0.:
            logger.writeln("AMBER prior: NOTE: chirality, NCS/stacking restraints and VDW repulsion against "
                           "symmetry copies are scaled by the same factor and have no AMBER equivalent")
        return int(serials.size)

    def calc_target_and_grad(self, target_only=False):
        self._context.setPositions(self._positions_nm())
        state = self._context.getState(getEnergy=True, getForces=(not target_only))
        energy = state.getPotentialEnergy().value_in_unit(self._unit.kilojoule_per_mole)
        if target_only:
            return float(energy), None

        forces = state.getForces(asNumpy=True).value_in_unit(self._unit.kilojoule_per_mole / self._unit.nanometer)
        # grad = dE/dx(A) = -force(nm) * 0.1
        grad_xyz = -numpy.array(forces) * 0.1
        grad_param = xyz_grad_to_param_grad(grad_xyz, self._atom_to_param_x, self.refine_params.n_params())
        return float(energy), grad_param
