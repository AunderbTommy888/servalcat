"""
Force-field templates for residues that the AMBER protein/nucleic-acid/water
force fields do not cover (ligands, ions with parameters, modified residues that
are not covalently linked).

The chemistry (elements, bond orders, formal charges, hydrogens) is taken from
the monomer library entry that Servalcat already uses for restraints, converted
to an RDKit molecule and then to an OpenFF Molecule. Parameters are assigned by
openmmforcefields (SMIRNOFF/OpenFF or GAFF template generators).

Optional dependencies: rdkit, openff-toolkit, openmmforcefields (and AmberTools
for AM1-BCC charges or GAFF; openff-nagl for NAGL charges).
"""
from __future__ import absolute_import, division, print_function, generators
import gemmi
from servalcat.utils import logger

CHARGE_METHODS = ("nagl", "am1bcc", "gasteiger")
DEFAULT_CHARGE_METHOD = "nagl"
DEFAULT_LIGAND_FF = "openff-2.2.1"


def _import_openff():
    try:
        from rdkit import Chem  # noqa: F401
        from openff.toolkit import Molecule  # noqa: F401
        import openmmforcefields  # noqa: F401
    except ImportError as e:
        raise RuntimeError(
            "Ligand force-field assignment needs rdkit, openff-toolkit and openmmforcefields "
            "(conda-forge packages; see docs/dev/build_uv_memo_ja.md). Original error: {}".format(e)) from e


def chemcomp_to_rdkit(cc, atom_names=None):
    """Build an RDKit molecule from a gemmi ChemComp.

    atom_names: restrict to these atom names (order preserved); default all atoms.
    Bond orders and formal charges are taken from the dictionary. Aromatic bonds are
    expected in Kekule form (as written by AceDRG); Deloc/Unspec bonds are first tried as
    single bonds and, if the valence model fails, bond orders are re-derived by RDKit's
    DetermineBondOrders from connectivity and total charge.
    Returns (mol, names) where names[i] is the atom name of RDKit atom i.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    names = [a.id for a in cc.atoms] if atom_names is None else list(atom_names)
    index = {n: i for i, n in enumerate(names)}
    cc_atoms = {a.id: a for a in cc.atoms}
    missing = [n for n in names if n not in cc_atoms]
    if missing:
        raise RuntimeError("atoms {} not found in monomer library entry {}".format(missing, cc.name))

    order_map = {gemmi.BondType.Single: Chem.BondType.SINGLE,
                 gemmi.BondType.Double: Chem.BondType.DOUBLE,
                 gemmi.BondType.Triple: Chem.BondType.TRIPLE,
                 gemmi.BondType.Aromatic: Chem.BondType.AROMATIC}

    def build(force_single):
        rw = Chem.RWMol()
        for n in names:
            a = cc_atoms[n]
            ra = Chem.Atom(a.el.atomic_number)
            ra.SetFormalCharge(int(round(a.charge)))
            ra.SetNoImplicit(True)
            rw.AddAtom(ra)
        undetermined = False
        for b in cc.rt.bonds:
            i, j = index.get(b.id1.atom), index.get(b.id2.atom)
            if i is None or j is None:
                continue
            if b.type == gemmi.BondType.Metal:
                continue  # coordination, not a covalent bond
            bt = order_map.get(b.type)
            if bt is None or force_single:
                bt = Chem.BondType.SINGLE
                if b.type not in order_map:
                    undetermined = True
            rw.AddBond(i, j, bt)
            if bt == Chem.BondType.AROMATIC:
                rw.GetAtomWithIdx(i).SetIsAromatic(True)
                rw.GetAtomWithIdx(j).SetIsAromatic(True)
                rw.GetBondBetweenAtoms(i, j).SetIsAromatic(True)
        return rw.GetMol(), undetermined

    mol, undetermined = build(force_single=False)
    try:
        if undetermined:
            raise ValueError("dictionary has Deloc/Unspec bonds")
        Chem.SanitizeMol(mol)
    except Exception as e1:
        # fall back: connectivity only + total charge -> bond orders
        logger.writeln(" ligand {}: re-deriving bond orders from connectivity ({})".format(cc.name, e1))
        mol, _ = build(force_single=True)
        total_charge = int(round(sum(cc_atoms[n].charge for n in names)))
        try:
            rdDetermineBonds.DetermineBondOrders(mol, charge=total_charge)
            Chem.SanitizeMol(mol)
        except Exception as e2:
            raise RuntimeError("Could not build a valid molecule for {} from the monomer library "
                               "(try --amber_ligand_smiles {}=<SMILES>): {}".format(cc.name, cc.name, e2)) from e2
    return mol, names


def chemcomp_to_openff_molecule(cc, atom_names=None, smiles=None):
    """OpenFF Molecule for a monomer library entry.

    If smiles is given it defines the chemistry (must include all hydrogens present in
    the model, i.e. same number of atoms per element as the residue)."""
    from openff.toolkit import Molecule
    from rdkit import Chem
    if smiles is not None:
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True, hydrogens_are_explicit=False)
        mol.name = cc.name
        return mol
    rdmol, names = chemcomp_to_rdkit(cc, atom_names)
    for i, n in enumerate(names):
        rdmol.GetAtomWithIdx(i).SetProp("_Name", n)
    # stereo from dictionary coordinates when available; otherwise leave undefined
    has_xyz = cc.has_coordinates() if callable(cc.has_coordinates) else cc.has_coordinates
    if has_xyz:
        conf = Chem.Conformer(rdmol.GetNumAtoms())
        for i, n in enumerate(names):
            p = cc.get_atom(n).xyz
            conf.SetAtomPosition(i, (p.x, p.y, p.z))
        rdmol.AddConformer(conf, assignId=True)
        Chem.AssignStereochemistryFrom3D(rdmol)
    mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    mol.name = cc.name
    return mol


def assign_ligand_charges(mol, method):
    """Partial charges for an OpenFF Molecule."""
    if method not in CHARGE_METHODS:
        raise RuntimeError("Unknown ligand charge method: {} (choose from {})".format(method, CHARGE_METHODS))
    if method == "nagl":
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
            from openff.nagl_models import list_available_nagl_models
        except ImportError as e:
            raise RuntimeError("--amber_ligand_charge nagl needs openff-nagl and openff-nagl-models") from e
        models = [m for m in list_available_nagl_models() if "am1bcc" in str(m)]
        if not models:
            raise RuntimeError("no NAGL am1bcc model found")
        model = sorted(models)[-1]
        mol.assign_partial_charges(str(model), toolkit_registry=NAGLToolkitWrapper())
        return
    if method == "gasteiger":
        from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
        mol.assign_partial_charges("gasteiger", toolkit_registry=RDKitToolkitWrapper())
        return
    # am1bcc: AmberTools (sqm must be on PATH) via openff-toolkit, or OpenEye if licensed
    try:
        mol.assign_partial_charges("am1bcc")
    except Exception as e:
        raise RuntimeError("AM1-BCC charge assignment failed (AmberTools 'sqm' must be on PATH; "
                           "alternatives: --amber_ligand_charge nagl or gasteiger). Original error: {}".format(e)) from e


def add_residue_bonds(topology, residue, cc):
    """Add intra-residue bonds of a non-standard residue to an OpenMM topology from the
    monomer library (PDBFile only knows bonds of standard residues)."""
    atoms = {a.name: a for a in residue.atoms()}
    existing = set()
    for a1, a2 in topology.bonds():
        if a1.residue is residue and a2.residue is residue:
            existing.add(frozenset((a1.index, a2.index)))
    n_added = 0
    for b in cc.rt.bonds:
        if b.type == gemmi.BondType.Metal:
            continue
        a1, a2 = atoms.get(b.id1.atom), atoms.get(b.id2.atom)
        if a1 is None or a2 is None:
            continue
        key = frozenset((a1.index, a2.index))
        if key in existing:
            continue
        topology.addBond(a1, a2)
        existing.add(key)
        n_added += 1
    return n_added


def find_ligand_residues(forcefield, topology, standard_residues=None):
    """Residues that need a generated template: non-standard residue names for which the
    ForceField has no template of the same name (ions, caps etc. are covered by the AMBER XMLs).
    ForceField.getUnmatchedResidues() is not used because it cannot ignore external bonds
    (a C-terminal residue without OXT would be reported)."""
    from openmm.app.pdbfile import PDBFile
    if standard_residues is None:
        standard_residues = set(PDBFile._standardResidues)
    templates = getattr(forcefield, "_templates", {})
    return [res for res in topology.residues()
            if res.name not in standard_residues and res.name not in templates]


def register_ligand_templates(forcefield, topology, monlib, ligand_ff=DEFAULT_LIGAND_FF,
                              charge_method="nagl", smiles_overrides=None, standard_residues=None):
    """Build OpenFF molecules for residues without AMBER templates from the monomer library
    and register a template generator on the OpenMM ForceField.

    Returns list of (residue name, n_atoms, total formal charge, ligand_ff) for registered molecules.
    Raises RuntimeError with a specific message when a residue cannot be handled.
    """
    from openmm.app.pdbfile import PDBFile
    smiles_overrides = smiles_overrides or {}
    if standard_residues is None:
        standard_residues = set(PDBFile._standardResidues)

    unmatched = find_ligand_residues(forcefield, topology, standard_residues)
    if not unmatched:
        return []
    _import_openff()
    # complete the connectivity of these residues first (needed for graph matching)
    for res in unmatched:
        cc = monlib.monomers[res.name] if res.name in monlib.monomers else None
        if cc is None:
            raise RuntimeError("residue {} {} is not in the monomer library; cannot build a force-field template".format(
                res.name, res.id))
        add_residue_bonds(topology, res, cc)

    # inter-residue bonds involving a non-standard residue are not supported (covalent ligands, modified residues)
    for a1, a2 in topology.bonds():
        if a1.residue is not a2.residue and (a1.residue.name not in standard_residues or a2.residue.name not in standard_residues):
            raise RuntimeError("residue {} {} is covalently linked to another residue; linked non-standard residues "
                               "are not supported by the ligand force-field assignment yet".format(
                                   a1.residue.name, a1.residue.id))

    names = []
    for res in unmatched:
        if res.name not in names:
            names.append(res.name)

    molecules = []
    info = []
    for name in names:
        cc = monlib.monomers[name]
        res = next(r for r in topology.residues() if r.name == name)
        atom_names = [a.name for a in res.atoms()]
        mol = chemcomp_to_openff_molecule(cc, atom_names, smiles=smiles_overrides.get(name))
        if mol.n_atoms != len(atom_names):
            raise RuntimeError("ligand {}: molecule has {} atoms but residue has {} (protonation/SMILES mismatch)".format(
                name, mol.n_atoms, len(atom_names)))
        logger.writeln(" ligand {}: {} atoms, formal charge {:+d}, assigning {} charges".format(
            name, mol.n_atoms, int(round(mol.total_charge.m)), charge_method))
        assign_ligand_charges(mol, charge_method)
        molecules.append(mol)
        info.append((name, mol.n_atoms, int(round(mol.total_charge.m)), ligand_ff))

    if ligand_ff.lower().startswith("gaff"):
        from openmmforcefields.generators import GAFFTemplateGenerator
        gen = GAFFTemplateGenerator(molecules=molecules, forcefield=ligand_ff)
    else:
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        gen = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield=ligand_ff)
    forcefield.registerTemplateGenerator(gen.generator)
    return info


def parse_smiles_overrides(items):
    """['BTN=C(...)', ...] -> {'BTN': 'C(...)'}"""
    out = {}
    for it in items or []:
        if "=" not in it:
            raise RuntimeError("--amber_ligand_smiles expects NAME=SMILES, got {}".format(it))
        k, v = it.split("=", 1)
        out[k.strip()] = v.strip()
    return out
