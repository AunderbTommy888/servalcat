from __future__ import absolute_import, division, print_function, generators
import unittest
import numpy
import os

from servalcat.refine import refine_spa
from servalcat.refine.ff_amber import xyz_grad_to_param_grad, AmberFFPrior
from servalcat.refine.refine import RefineParams
from servalcat import utils
from servalcat import ext


class TestAmberPhase1(unittest.TestCase):
    @staticmethod
    def _renumber_serials(st):
        for i, cra in enumerate(st[0].all()):
            cra.atom.serial = i + 1

    @staticmethod
    def _remove_waters(st):
        for chain in st[0]:
            for i in reversed(range(len(chain))):
                if chain[i].name == "HOH":
                    del chain[i]

    def test_parse_args_amber_defaults(self):
        args = refine_spa.parse_args([
            "--halfmaps", "half1.mrc", "half2.mrc",
            "--model", "model.pdb",
            "-d", "3.0",
        ])
        self.assertFalse(args.amber_enable)
        self.assertAlmostEqual(args.amber_weight, 0.1)
        self.assertEqual(args.amber_nonbonded, "NoCutoff")
        self.assertEqual(args.amber_platform, "Reference")
        self.assertEqual(args.amber_hessian_diag, 1000.0)
        self.assertEqual(args.amber_hessian_mode, "bonded")
        self.assertEqual(args.amber_his_state, "HIP")
        self.assertEqual(args.amber_forcefield, ["amber14-all.xml", "amber14/tip3p.xml"])

    def test_parse_args_amber_override(self):
        args = refine_spa.parse_args([
            "--halfmaps", "half1.mrc", "half2.mrc",
            "--model", "model.pdb",
            "-d", "3.0",
            "--amber_enable",
            "--amber_weight", "0.25",
            "--amber_nonbonded", "CutoffNonPeriodic",
            "--amber_platform", "CPU",
            "--amber_hessian_diag", "7.5",
            "--amber_hessian_mode", "const",
            "--amber_his_state", "HIE",
            "--amber_forcefield", "amber14-all.xml", "amber14/tip3p.xml",
        ])
        self.assertTrue(args.amber_enable)
        self.assertAlmostEqual(args.amber_weight, 0.25)
        self.assertEqual(args.amber_nonbonded, "CutoffNonPeriodic")
        self.assertEqual(args.amber_platform, "CPU")
        self.assertAlmostEqual(args.amber_hessian_diag, 7.5)
        self.assertEqual(args.amber_hessian_mode, "const")
        self.assertEqual(args.amber_his_state, "HIE")

    def test_xyz_grad_to_param_grad(self):
        grad_xyz = numpy.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        # atom 1 is excluded from xyz refinement
        atom_to_param = [0, -1, 2]
        g = xyz_grad_to_param_grad(grad_xyz, atom_to_param, n_params=9)
        expected = numpy.array([
            1.0, 2.0, 3.0,
            0.0, 0.0, 0.0,
            7.0, 8.0, 9.0,
        ])
        numpy.testing.assert_allclose(g, expected)

    def test_amber_prior_energy_grad_optional(self):
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")

        st = utils.fileio.read_structure("tests/1l2h/1l2h.cif.gz")
        self._remove_waters(st)
        self._renumber_serials(st)
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set; needed to add hydrogens for OpenMM templates")
        monlib = utils.restraints.load_monomer_library(st)
        utils.restraints.add_hydrogens(st, monlib, "elec")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        prior = AmberFFPrior(st, rp,
                             forcefield_files=["amber14-all.xml", "amber14/tip3p.xml"],
                             nonbonded_method="NoCutoff",
                             platform_name="Reference",
                             hessian_diag=10.0)
        # 1l2h contains alternative conformations; only the first conformer goes to OpenMM
        n_total = sum(1 for _ in st[0].all())
        self.assertGreater(prior.n_excluded_atoms, 0)
        self.assertEqual(prior.n_ff_atoms + prior.n_excluded_atoms, n_total)
        ene, grad = prior.calc_target_and_grad(target_only=False)
        self.assertTrue(numpy.isfinite(ene))
        self.assertEqual(len(grad), rp.n_params())
        self.assertTrue(numpy.all(numpy.isfinite(grad)))
        # excluded altloc atoms must receive zero force-field gradient and zero Hessian diagonal
        hdiag = prior.hessian_diag_vector()
        self.assertEqual(len(hdiag), rp.n_params())
        atom_to_param = list(rp.atom_to_param(ext.RefineParams.Type.X))
        for cra in st[0].all():
            pidx = atom_to_param[cra.atom.serial - 1]
            self.assertGreaterEqual(pidx, 0)
            if cra.atom.altloc not in ('\0', 'A'):
                numpy.testing.assert_array_equal(grad[pidx*3:pidx*3+3], 0.0)
                numpy.testing.assert_array_equal(hdiag[pidx*3:pidx*3+3], 0.0)
            else:
                # bonded estimate is never below the floor and is of the expected magnitude (kJ/mol/A^2)
                self.assertTrue(numpy.all(hdiag[pidx*3:pidx*3+3] >= prior.hessian_diag))
                self.assertTrue(numpy.all(hdiag[pidx*3:pidx*3+3] < 5e4))
        # const mode reproduces the old behaviour (same value for every force-field atom)
        prior_c = AmberFFPrior(st, rp, platform_name="Reference", hessian_diag=10.0, hessian_mode="const")
        hc = prior_c.hessian_diag_vector()
        self.assertEqual(set(numpy.unique(hc[hc > 0]).tolist()), {10.0})

    def test_bonded_hessian_diag_matches_finite_difference(self):
        # water-like 3-atom system at equilibrium: Gauss-Newton diagonal equals the exact Hessian diagonal
        from servalcat.refine.ff_amber import bonded_hessian_diag
        r0, t0, kb, kt = 0.9572, numpy.deg2rad(104.52), 4627.5, 836.8
        pos = numpy.array([[0., 0., 0.],
                           [r0, 0., 0.],
                           [r0 * numpy.cos(t0), r0 * numpy.sin(t0), 0.]])
        bonds = (numpy.array([1, 2]), numpy.array([0, 0]), numpy.array([kb, kb]))
        angles = (numpy.array([1]), numpy.array([0]), numpy.array([2]), numpy.array([kt]))
        def energy(p):
            e = 0.
            for i, j, k in zip(*bonds):
                e += 0.5 * k * (numpy.linalg.norm(p[i] - p[j]) - r0)**2
            for i, j, k, kth in zip(*angles):
                a, b = p[i] - p[j], p[k] - p[j]
                t = numpy.arccos(numpy.dot(a, b) / numpy.linalg.norm(a) / numpy.linalg.norm(b))
                e += 0.5 * kth * (t - t0)**2
            return e
        est = bonded_hessian_diag(pos, bonds, angles)
        h = 1e-4
        for i in range(3):
            for ax in range(3):
                pp, pm = pos.copy(), pos.copy()
                pp[i, ax] += h; pm[i, ax] -= h
                fd = (energy(pp) - 2 * energy(pos) + energy(pm)) / h**2
                self.assertAlmostEqual(est[i, ax], fd, delta=max(1.0, 1e-3 * abs(fd)))

    def test_his_state_excludes_protons(self):
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        st = utils.fileio.read_structure("tests/1l2h/1l2h.cif.gz")
        self._remove_waters(st)
        self._renumber_serials(st)
        monlib = utils.restraints.load_monomer_library(st)
        utils.restraints.add_hydrogens(st, monlib, "elec")
        self._renumber_serials(st)
        n_hd1 = sum(1 for cra in st[0].all() if cra.residue.name == "HIS" and cra.atom.name == "HD1")
        n_he2 = sum(1 for cra in st[0].all() if cra.residue.name == "HIS" and cra.atom.name == "HE2")
        self.assertGreater(n_hd1, 0)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        base = AmberFFPrior(st, rp, platform_name="Reference", his_state="HIP")
        hie = AmberFFPrior(st, rp, platform_name="Reference", his_state="HIE")
        hid = AmberFFPrior(st, rp, platform_name="Reference", his_state="HID")
        self.assertEqual(hie.n_ff_atoms, base.n_ff_atoms - n_hd1)
        self.assertEqual(hid.n_ff_atoms, base.n_ff_atoms - n_he2)
        # HIS residues get the HIE/HID templates: all three systems build and give finite energies
        for p in (base, hie, hid):
            ene, grad = p.calc_target_and_grad()
            self.assertTrue(numpy.isfinite(ene))
            self.assertTrue(numpy.all(numpy.isfinite(grad)))
        with self.assertRaises(RuntimeError):
            AmberFFPrior(st, rp, platform_name="Reference", his_state="HIX")

    def test_amber_prior_unsupported_residue_error(self):
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")

        st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        with self.assertRaises(RuntimeError) as cm:
            AmberFFPrior(st, rp,
                         forcefield_files=["amber14-all.xml", "amber14/tip3p.xml"],
                         nonbonded_method="NoCutoff",
                         platform_name="Reference",
                         hessian_diag=10.0)
        self.assertIn("OpenMM failed to build an AMBER system", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
