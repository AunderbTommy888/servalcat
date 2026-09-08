from __future__ import absolute_import, division, print_function, generators
import unittest
import numpy
import os

from servalcat.refine import refine_spa
from servalcat.refine.ff_amber import xyz_grad_to_param_grad, AmberFFPrior
from servalcat.refine.refine import RefineParams
from servalcat import utils


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
        self.assertEqual(args.amber_hessian_diag, 10.0)
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
            "--amber_forcefield", "amber14-all.xml", "amber14/tip3p.xml",
        ])
        self.assertTrue(args.amber_enable)
        self.assertAlmostEqual(args.amber_weight, 0.25)
        self.assertEqual(args.amber_nonbonded, "CutoffNonPeriodic")
        self.assertEqual(args.amber_platform, "CPU")
        self.assertAlmostEqual(args.amber_hessian_diag, 7.5)

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
        ene, grad = prior.calc_target_and_grad(target_only=False)
        self.assertTrue(numpy.isfinite(ene))
        self.assertEqual(len(grad), rp.n_params())
        self.assertTrue(numpy.all(numpy.isfinite(grad)))

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
