"""
AMBER force-field prior for Servalcat refinement (Phase 1).

This module is intentionally minimal:
- Optional dependency on OpenMM
- Energy + Cartesian gradient contribution only
- XYZ refinement only (no ADP/occ coupling)
"""
from __future__ import absolute_import, division, print_function, generators
import os
import tempfile
import numpy
import gemmi
from servalcat import ext


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


class AmberFFPrior:
    def __init__(self, st, refine_params,
                 forcefield_files=None,
                 nonbonded_method="NoCutoff",
                 platform_name="Reference",
                 hessian_diag=10.):
        self.st = st
        self.refine_params = refine_params
        self.hessian_diag = float(hessian_diag)
        self._openmm = None
        self._unit = None
        self._context = None
        self._integrator = None
        self._n_atoms = len(self.refine_params.atoms)
        self._ordered_atoms = [cra.atom for cra in self.st[0].all()]

        if forcefield_files is None:
            forcefield_files = ["amber14-all.xml", "amber14/tip3p.xml"]
        self.forcefield_files = list(forcefield_files)
        self.nonbonded_method = nonbonded_method
        self.platform_name = platform_name

        self._setup_openmm_system()

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

        with tempfile.TemporaryDirectory(prefix="servalcat_amber_") as dtmp:
            pdb_path = os.path.join(dtmp, "amber_input.pdb")
            self.st.write_pdb(pdb_path)
            pdb = app.PDBFile(pdb_path)

            ff = app.ForceField(*self.forcefield_files)
            if not hasattr(app, self.nonbonded_method):
                raise RuntimeError("Unknown OpenMM nonbonded method: {}".format(self.nonbonded_method))
            nb_method = getattr(app, self.nonbonded_method)
            try:
                system = ff.createSystem(pdb.topology,
                                         nonbondedMethod=nb_method,
                                         constraints=None,
                                         rigidWater=False,
                                         ignoreExternalBonds=True)
            except Exception as e:
                raise RuntimeError(
                    "OpenMM failed to build an AMBER system. "
                    "This usually means the model contains residues without force-field templates "
                    "(e.g. ligands) and needs additional parameter XML files.\n"
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

        if self._topology_n_atoms != len(self._ordered_atoms):
            raise RuntimeError(
                "OpenMM topology atom count ({}) does not match structure atom count ({})".format(
                    self._topology_n_atoms, len(self._ordered_atoms)
                )
            )

        atom_to_param_ref = list(self.refine_params.atom_to_param(ext.RefineParams.Type.X))
        self._atom_to_param_x = []
        for i, _ in enumerate(self._topology_atoms):
            pidx = atom_to_param_ref[i] if i < len(atom_to_param_ref) else -1
            self._atom_to_param_x.append(pidx)

    def _positions_nm(self):
        # OpenMM expects nm; gemmi positions are angstrom.
        out = []
        for atom in self._ordered_atoms:
            p = atom.pos
            out.append(self._openmm.Vec3(p.x * 0.1, p.y * 0.1, p.z * 0.1))
        return out

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
