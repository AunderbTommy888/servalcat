from __future__ import absolute_import, division, print_function, generators
import unittest
import numpy
import os
import sys

from servalcat.refine import refine_spa
from servalcat.refine.ff_amber import xyz_grad_to_param_grad, AmberFFPrior
from servalcat.refine.refine import RefineParams
from servalcat import utils
from servalcat import ext
import gemmi


class FakePrior:
    """Minimal stand-in for AmberFFPrior: a quadratic well in parameter space, read through the
    RefineParams object exactly like the real prior reads the gemmi atoms. No OpenMM needed."""
    hessian_diag = 1000.
    hessian_mode = "const"
    his_state = "HIP"

    def __init__(self, params, k=1000., offset=0.1):
        self.params = params
        self.n = params.n_params()
        self.k = k
        # the minimum sits away from the current parameters, so the gradient is non-zero there
        self.x0 = numpy.array(params.get_x()) - offset

    def _d(self):
        return numpy.array(self.params.get_x()) - self.x0

    def calc_target_and_grad(self, target_only=False):
        d = self._d()
        f = 0.5 * self.k * float(numpy.dot(d, d))
        return f, (None if target_only else self.k * d)

    def hessian_diag_vector(self):
        return numpy.full(self.n, self.k)

    def hessian_matrix(self):
        import scipy.sparse
        off = 0.1 * self.k  # small symmetric coupling between neighbouring parameters
        m = scipy.sparse.diags([numpy.full(self.n - 1, off), numpy.full(self.n, float(self.k)),
                                numpy.full(self.n - 1, off)], [-1, 0, 1], format="csr")
        return m


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
        self.assertEqual(args.amber_nonbonded, "NoCutoff")
        self.assertEqual(args.amber_platform, "Reference")
        self.assertEqual(args.amber_hessian_diag, 1000.0)
        self.assertEqual(args.amber_replace_geom, 1.0)   # full replacement by default
        self.assertIsNone(args.amber_weight)
        self.assertEqual(refine_spa.resolve_amber_weight(args), 1.0)
        self.assertEqual(refine_spa.resolve_amber_weight(
            refine_spa.parse_args(["--halfmaps", "a", "b", "--model", "m.pdb", "-d", "3.0",
                                   "--amber_replace_geom", "0"])), 0.1)
        self.assertIsNone(args.amber_weight_auto)
        self.assertEqual(refine_spa.parse_args(["--halfmaps", "a", "b", "--model", "m.pdb", "-d", "3.0",
                                                "--amber_weight_auto"]).amber_weight_auto, 0.3)
        self.assertEqual(refine_spa.parse_args(["--halfmaps", "a", "b", "--model", "m.pdb", "-d", "3.0",
                                                "--amber_weight_auto", "0.5"]).amber_weight_auto, 0.5)
        self.assertEqual(args.amber_hessian_mode, "bonded")
        self.assertFalse(args.amber_hessian_offdiag)
        self.assertEqual(args.amber_minimizer, "gn")
        self.assertEqual(args.amber_lbfgs_maxiter, 20)
        self.assertEqual(args.amber_lm_damping, 0.)
        self.assertIsNone(args.amber_ls_trials)
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
            "--amber_minimizer", "lbfgs",
            "--amber_lbfgs_maxiter", "7",
            "--amber_his_state", "HIE",
            "--amber_forcefield", "amber14-all.xml", "amber14/tip3p.xml",
        ])
        self.assertTrue(args.amber_enable)
        self.assertAlmostEqual(args.amber_weight, 0.25)
        self.assertAlmostEqual(refine_spa.resolve_amber_weight(args), 0.25)
        self.assertEqual(args.amber_nonbonded, "CutoffNonPeriodic")
        self.assertEqual(args.amber_platform, "CPU")
        self.assertAlmostEqual(args.amber_hessian_diag, 7.5)
        self.assertEqual(args.amber_hessian_mode, "const")
        self.assertEqual(args.amber_minimizer, "lbfgs")
        self.assertEqual(args.amber_lbfgs_maxiter, 7)
        self.assertEqual(args.amber_his_state, "HIE")

    def test_check_amber_args(self):
        base = ["--halfmaps", "half1.mrc", "half2.mrc", "--model", "model.pdb", "-d", "3.0", "--amber_enable"]
        refine_spa.check_amber_args(refine_spa.parse_args(base))  # --hydrogen all is the default
        for extra in (["--hydrogen", "no"], ["--hydrogen", "yes"], ["--unrestrained"], ["--jellyonly"]):
            with self.assertRaises(SystemExit):
                refine_spa.check_amber_args(refine_spa.parse_args(base + extra))
        # without --amber_enable nothing is checked
        refine_spa.check_amber_args(refine_spa.parse_args(base[:-1] + ["--hydrogen", "no"]))
        with self.assertRaises(SystemExit):  # PME is no longer a valid choice
            refine_spa.parse_args(base + ["--amber_nonbonded", "PME"])
        # the replacement ratio must be a fraction, and it is incompatible with the automatic weight
        for extra in (["--amber_replace_geom", "1.5"], ["--amber_replace_geom", "-0.1"],
                      ["--amber_weight_auto", "0.3"], ["--unrestrained"]):
            with self.assertRaises(SystemExit, msg=str(extra)):
                refine_spa.check_amber_args(refine_spa.parse_args(base + extra))
        # with the ratio at 0 the automatic weight is allowed again
        refine_spa.check_amber_args(refine_spa.parse_args(base + ["--amber_replace_geom", "0",
                                                                  "--amber_weight_auto", "0.3"]))
        # off-diagonal Hessian needs the bonded mode
        refine_spa.check_amber_args(refine_spa.parse_args(base + ["--amber_hessian_offdiag"]))
        with self.assertRaises(SystemExit):
            refine_spa.check_amber_args(refine_spa.parse_args(
                base + ["--amber_hessian_offdiag", "--amber_hessian_mode", "const"]))
        # minimiser / off-diagonal options require --amber_enable
        for extra in (["--amber_minimizer", "lbfgs"], ["--amber_hessian_offdiag"]):
            with self.assertRaises(SystemExit):
                refine_spa.check_amber_args(refine_spa.parse_args(base[:-1] + extra))

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

    def test_determine_ff_weight(self):
        # geometry-only Refine object with a fake force-field prior: the weight must scale the
        # ff gradient norm to RATIO x geometry gradient norm over xyz parameters (no OpenMM needed)
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine.refine import Refine, Geom, load_config
        st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
        monlib = utils.restraints.load_monomer_library(st)
        topo, _ = utils.restraints.prepare_topology(st, monlib, h_change=gemmi.HydrogenChange.NoChange)
        self._renumber_serials(st)
        from servalcat.refmac import refmac_keywords
        params = refmac_keywords.parse_keywords([])
        cfg = load_config(None, None, params)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=cfg)
        geom = Geom(st, topo, monlib, rp, cfg, shake_rms=0.3, params=params)
        geom.setup_nonbonded()
        geom.setup_target()

        prior = FakePrior(rp)
        ref = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=0.1, ff_weight_auto=0.5)
        self.assertFalse(ref.ff_weight_determined)
        ref.calc_target(1)
        self.assertTrue(ref.ff_weight_determined)
        sel = rp.vec_selection(ext.RefineParams.Type.X)
        g_geom = numpy.linalg.norm(numpy.array(geom.geom.target.vn)[sel])
        g_ff = numpy.linalg.norm(prior.calc_target_and_grad()[1][sel])
        self.assertAlmostEqual(ref.ff_weight * g_ff / g_geom, 0.5, places=6)
        w1 = ref.ff_weight
        ref.calc_target(1)  # determined once only
        self.assertEqual(ref.ff_weight, w1)

    def _geom_setup(self, shake=0.3):
        from servalcat.refine.refine import Geom, load_config
        from servalcat.refmac import refmac_keywords
        st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
        monlib = utils.restraints.load_monomer_library(st)
        topo, _ = utils.restraints.prepare_topology(st, monlib, h_change=gemmi.HydrogenChange.NoChange)
        self._renumber_serials(st)
        params = refmac_keywords.parse_keywords([])
        cfg = load_config(None, None, params)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=cfg)
        geom = Geom(st, topo, monlib, rp, cfg, shake_rms=shake, params=params)
        geom.setup_nonbonded()
        geom.setup_target()
        return st, rp, geom, cfg

    def test_run_cycle_lbfgs(self):
        """L-BFGS version: diagonal-only Hessian used as preconditioner, curvature from gradient history."""
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine.refine import Refine
        st, rp, geom, cfg = self._geom_setup()
        prior = FakePrior(rp)
        ref = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=0.05,
                     minimizer="lbfgs", lbfgs_maxiter=15)
        x0 = ref.get_x()
        f0 = ref.calc_target(1)
        ok, scale, f1 = ref.run_cycle(weight=1)
        self.assertTrue(ok)
        self.assertLess(f1, f0)
        self.assertGreater(numpy.abs(ref.get_x() - x0).max(), 1e-6)
        # the trust region of scale_shifts must be respected
        self.assertLessEqual(numpy.abs(ref.get_x() - x0).max(), 1.0 + 1e-9)
        self.assertFalse(ref.quiet)  # restored after the inner iterations

    def test_lm_damping_makes_matrix_positive_definite(self):
        """A positive Levenberg-Marquardt ridge lifts the rigid-body null modes, so the normal
        matrix is strictly positive definite instead of only semi-definite."""
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine.refine import Refine
        st, rp, geom, cfg = self._geom_setup(shake=0.0)
        prior = FakePrior(rp)
        geom.calc_target(False)
        Hg = geom.geom.target.am_spmat.toarray()
        ev0 = numpy.linalg.eigvalsh((Hg + Hg.T) / 2)
        self.assertGreater(ev0.min(), -1e-6 * abs(ev0).max())   # restraints alone: PSD
        for lam in (10., 100.):
            A = Hg + lam * numpy.eye(rp.n_params())
            self.assertGreater(numpy.linalg.eigvalsh((A + A.T) / 2).min(), lam * 0.99)
        # a negative ridge is refused
        with self.assertRaises(RuntimeError):
            Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=0.05, lm_damping=-1.)
        # and a cycle with the ridge still decreases the target
        ref = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=0.05,
                     ff_offdiag=True, lm_damping=100.)
        f0 = ref.calc_target(1)
        ok, _, f1 = ref.run_cycle(weight=1)
        self.assertTrue(ok)
        self.assertLess(f1, f0)

    def test_run_cycle_offdiag_hessian(self):
        """Off-diagonal version: same Gauss-Newton cycle, full force-field Hessian in the matrix."""
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine.refine import Refine
        res = {}
        for offdiag in (False, True):
            st, rp, geom, cfg = self._geom_setup()
            prior = FakePrior(rp)
            ref = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=0.05, ff_offdiag=offdiag)
            x0 = ref.get_x()
            f0 = ref.calc_target(1)
            ok, scale, f1 = ref.run_cycle(weight=1)
            self.assertTrue(ok)
            self.assertLess(f1, f0)
            res[offdiag] = ref.get_x() - x0
        # the off-diagonal terms must change the step
        self.assertGreater(numpy.abs(res[True] - res[False]).max(), 1e-6)

    def test_bonded_hessian_is_positive_semidefinite(self):
        """The assembled bonded Hessian must never have a negative eigenvalue: the normal matrix
        relies on it (adding a PSD matrix with a non-negative weight preserves positive definiteness).
        Synthetic geometries, no OpenMM needed."""
        from servalcat.refine.ff_amber import (assemble_bonded_hessian, bonded_hessian_diag,
                                               angle_gradients, MIN_SIN_ANGLE)
        rng = numpy.random.default_rng(0)
        cases = {}
        # water-like, bent
        r0, t0 = 0.9572, numpy.deg2rad(104.52)
        cases["bent"] = numpy.array([[0., 0., 0.], [r0, 0., 0.],
                                     [r0*numpy.cos(t0), r0*numpy.sin(t0), 0.]])
        # nearly linear (the pathological case for the 1/sin(theta) factor)
        for deg in (179.0, 179.999, 180.0):
            t = numpy.deg2rad(deg)
            cases["linear_%s" % deg] = numpy.array([[0., 0., 0.], [r0, 0., 0.],
                                                    [r0*numpy.cos(t), r0*numpy.sin(t), 0.]])
        # random cloud
        cases["random"] = rng.normal(size=(3, 3))
        bonds = (numpy.array([1, 2]), numpy.array([0, 0]), numpy.array([4627.5, 4627.5]))
        angles = (numpy.array([1]), numpy.array([0]), numpy.array([2]), numpy.array([836.8]))
        a2p = numpy.array([0, 1, 2])
        for name, pos in cases.items():
            A = assemble_bonded_hessian(pos, bonds, angles, a2p, 9).toarray()
            self.assertLess(numpy.abs(A - A.T).max(), 1e-8, name)
            ev = numpy.linalg.eigvalsh(A)
            self.assertGreater(ev.min(), -1e-6 * max(1., numpy.abs(A).max()), "%s: %s" % (name, ev[:3]))
            # the diagonal must agree with the per-atom estimate used by the diagonal-only version
            numpy.testing.assert_allclose(numpy.diag(A), bonded_hessian_diag(pos, bonds, angles).ravel(),
                                          rtol=1e-9, atol=1e-8)
            # rigid-body modes: translation and (for a rigid rotation) zero curvature
            for ax in range(3):
                t = numpy.zeros(9); t[ax::3] = 1.
                self.assertLess(abs(float(t @ A @ t)), 1e-6 * numpy.abs(A).max(), name)
            # a floor only adds a non-negative diagonal, so it cannot break PSD
            Af = assemble_bonded_hessian(pos, bonds, angles, a2p, 9, diag_floor=1000.).toarray()
            self.assertGreaterEqual(numpy.linalg.eigvalsh(Af).min(), -1e-6 * numpy.abs(Af).max())
            self.assertTrue(numpy.all(numpy.diag(Af) >= 1000. - 1e-6), name)
        # the 1/sin(theta) factor must be clamped, not infinite
        gi, gj, gk, n = angle_gradients(cases["linear_180.0"], numpy.array([1]), numpy.array([0]), numpy.array([2]))
        self.assertEqual(n, 1)
        self.assertTrue(numpy.all(numpy.isfinite(gi)) and numpy.all(numpy.isfinite(gk)))
        self.assertLessEqual(numpy.abs(gi).max(), 1. / (r0 * MIN_SIN_ANGLE) + 1e-9)

    def test_non_positive_force_constants_are_dropped(self):
        """A negative force constant would make the Gauss-Newton block negative semi-definite."""
        from servalcat.refine.ff_amber import _keep_positive
        bonds = (numpy.array([0, 1, 2]), numpy.array([1, 2, 3]), numpy.array([100., -5., 0.]))
        kept, n = _keep_positive(bonds, 2)
        self.assertEqual(n, 2)
        numpy.testing.assert_array_equal(kept[2], [100.])
        numpy.testing.assert_array_equal(kept[0], [0])
        angles = (numpy.array([0]), numpy.array([1]), numpy.array([2]), numpy.array([50.]))
        kept, n = _keep_positive(angles, 3)
        self.assertEqual(n, 0)
        self.assertIs(kept, angles)

    def test_negative_weights_rejected(self):
        base = ["--halfmaps", "a", "b", "--model", "m.pdb", "-d", "3.0", "--amber_enable"]
        for extra in (["--amber_weight", "-0.1"], ["--amber_weight_auto", "-1"],
                      ["--amber_hessian_diag", "-1"], ["--amber_ls_trials", "0"],
                      ["--amber_lm_damping", "-1"]):
            with self.assertRaises(SystemExit, msg=str(extra)):
                refine_spa.check_amber_args(refine_spa.parse_args(base + extra))
        with self.assertRaises(RuntimeError):  # also guarded in the prior itself
            from servalcat.refine.ff_amber import AmberFFPrior as P
            P.__init__.__wrapped__ if False else None
            st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
            self._renumber_serials(st)
            rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
            P(st, rp, platform_name="Reference", hessian_diag=-1.)

    def test_hessian_matrix_offdiag(self):
        """hessian_matrix(): symmetric, diagonal identical to hessian_diag_vector(), and the
        bonded Gauss-Newton Hessian keeps its rigid-body null mode when the floor is disabled."""
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        import scipy.sparse
        st = utils.fileio.read_structure("tests/1l2h/1l2h.cif.gz")
        self._remove_waters(st)
        self._renumber_serials(st)
        monlib = utils.restraints.load_monomer_library(st)
        utils.restraints.add_hydrogens(st, monlib, "nucl")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        prior = AmberFFPrior(st, rp, platform_name="Reference", his_state="HIE")
        A = prior.hessian_matrix()
        self.assertEqual(A.shape, (rp.n_params(), rp.n_params()))
        self.assertLess(abs(A - A.T).max(), 1e-6)
        numpy.testing.assert_allclose(A.diagonal(), prior.hessian_diag_vector(), rtol=1e-9, atol=1e-6)
        self.assertGreater(A.nnz, (A.diagonal() != 0).sum())  # off-diagonal entries exist
        # plain-loop reference
        B = scipy.sparse.lil_matrix(A.shape)
        for pa, pb, H in prior.hessian_blocks():
            same = numpy.array_equal(pa, pb)
            for m in range(len(pa)):
                i, j = int(pa[m]), int(pb[m])
                if i < 0 or j < 0:
                    continue
                B[i*3:i*3+3, j*3:j*3+3] += H[m]
                if not same:
                    B[j*3:j*3+3, i*3:i*3+3] += H[m].T
        B = B.tocsr()
        dv = prior.hessian_diag_vector()
        B = B + scipy.sparse.diags(numpy.maximum(dv - B.diagonal(), 0.), format="csr")
        self.assertLess(abs(A - B).max(), 1e-6)
        # rigid-body null mode with the floor disabled
        prior0 = AmberFFPrior(st, rp, platform_name="Reference", his_state="HIE", hessian_diag=0.)
        A0 = prior0.hessian_matrix()
        a2p = numpy.asarray(prior0._atom_to_param_x)
        for ax in range(3):
            t = numpy.zeros(rp.n_params())
            t[a2p[a2p >= 0]*3 + ax] = 1.
            self.assertLess(abs(float(t @ (A0 @ t))), 1e-5 * A0.diagonal().max())
        with self.assertRaises(RuntimeError):  # const mode has no off-diagonal Hessian
            AmberFFPrior(st, rp, platform_name="Reference", hessian_mode="const").hessian_matrix()

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

    def test_replace_geom_restraints(self):
        """The per-atom geometry weight is scaled for the atoms the force field covers, and atoms
        outside it keep their classical restraints."""
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        st = utils.fileio.read_structure("tests/1l2h/1l2h.cif.gz")   # has alternative conformers
        self._remove_waters(st)
        self._renumber_serials(st)
        monlib = utils.restraints.load_monomer_library(st)
        utils.restraints.add_hydrogens(st, monlib, "nucl")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        prior = AmberFFPrior(st, rp, platform_name="Reference", his_state="HIE")
        self.assertGreater(prior.n_excluded_atoms, 0)
        w = rp.geom_weights
        numpy.testing.assert_allclose(numpy.asarray(w), 1., rtol=0, atol=1e-6)
        self.assertIsNone(prior.geom_weights_full)
        n = prior.replace_geom_restraints(1.0)
        self.assertEqual(n, prior.n_ff_atoms)
        self.assertIsNotNone(prior.geom_weights_full)
        w = numpy.asarray(rp.geom_weights)
        ff_serials = numpy.array([a.serial for a in prior._ff_atoms]) - 1
        numpy.testing.assert_allclose(w[ff_serials], 0., atol=1e-7)
        excluded = numpy.setdiff1d(numpy.arange(len(w)), ff_serials)
        self.assertEqual(len(excluded), prior.n_excluded_atoms)
        numpy.testing.assert_allclose(w[excluded], 1., atol=1e-7)   # keep full restraints
        # a partial ratio scales instead of removing, and multiplies rather than overwriting
        rp2 = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        numpy.asarray(rp2.geom_weights)[:] = 0.5
        prior2 = AmberFFPrior(st, rp2, platform_name="Reference", his_state="HIE")
        prior2.replace_geom_restraints(0.4)
        w2 = numpy.asarray(rp2.geom_weights)
        numpy.testing.assert_allclose(w2[ff_serials], 0.5 * 0.6, rtol=1e-6)
        numpy.testing.assert_allclose(w2[excluded], 0.5, rtol=1e-6)
        for bad in (-0.1, 1.1):
            with self.assertRaises(RuntimeError):
                prior2.replace_geom_restraints(bad)
        self.assertEqual(prior2.replace_geom_restraints(0.), 0)   # a no-op

    def test_model_stats_restores_weights(self):
        """The statistics pass must report the classical geometry even when the restraints carry no
        weight in the target, and must leave the target weights untouched."""
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine.refine import Refine
        st, rp, geom, cfg = self._geom_setup(shake=0.3)
        prior = FakePrior(rp)
        full = numpy.array(rp.geom_weights, copy=True)
        numpy.asarray(rp.geom_weights)[:] = 0.
        ref = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=1.,
                     geom_weights_full=full)
        stats = ref.model_stats(show_outliers=False)
        self.assertIn("Bond distances, non H", stats["summary"]["r.m.s.Z"])
        numpy.testing.assert_allclose(numpy.asarray(rp.geom_weights), 0., atol=1e-7)
        # without the backup the zero-weight restraints drop out of the report
        ref2 = Refine(st, geom, cfg, rp, ll=None, ff_prior=prior, ff_weight=1.)
        self.assertNotIn("Bond distances, non H", ref2.model_stats(show_outliers=False)["summary"]["r.m.s.Z"])

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
        with self.assertRaises(RuntimeError):
            AmberFFPrior(st, rp, platform_name="Reference", nonbonded_method="PME")

    def test_amber_prior_unsupported_residue_error(self):
        try:
            import openmm  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM is not installed")

        st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        # ligand parameterisation disabled: OpenMM template error
        with self.assertRaises(RuntimeError) as cm:
            AmberFFPrior(st, rp,
                         forcefield_files=["amber14-all.xml", "amber14/tip3p.xml"],
                         nonbonded_method="NoCutoff",
                         platform_name="Reference",
                         hessian_diag=10.0, ligand_ff="none")
        self.assertIn("OpenMM failed to build an AMBER system", str(cm.exception))
        # ligand parameterisation enabled but no monomer library given
        with self.assertRaises(RuntimeError) as cm:
            AmberFFPrior(st, rp, platform_name="Reference")
        self.assertIn("no monomer library", str(cm.exception))

    def test_parse_args_ligand_options(self):
        args = refine_spa.parse_args(["--halfmaps", "a", "b", "--model", "m.pdb", "-d", "3.0"])
        self.assertEqual(args.amber_ligand_ff, "openff-2.2.1")
        self.assertEqual(args.amber_ligand_charge, "nagl")
        self.assertEqual(args.amber_ligand_smiles, [])
        from servalcat.refine import ff_ligand
        self.assertEqual(ff_ligand.parse_smiles_overrides(["BTN=C(=O)[O-]", "X=CC"]), {"BTN": "C(=O)[O-]", "X": "CC"})
        with self.assertRaises(RuntimeError):
            ff_ligand.parse_smiles_overrides(["BTN"])

    def test_chemcomp_to_rdkit(self):
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("rdkit is not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        from servalcat.refine import ff_ligand
        import gemmi
        ml = gemmi.read_monomer_lib(os.environ["CLIBD_MON"], ["BTN", "TRP", "SO4", "SPK"])
        for name, expected_charge in (("BTN", -1), ("SO4", -2), ("SPK", 4)):
            mol, names = ff_ligand.chemcomp_to_rdkit(ml.monomers[name])
            self.assertEqual(mol.GetNumAtoms(), len(ml.monomers[name].atoms))
            self.assertEqual(Chem.GetFormalCharge(mol), expected_charge)
            self.assertEqual(len(names), mol.GetNumAtoms())
        # aromatic ring (Kekule form in the dictionary) is perceived as aromatic by RDKit
        mol, names = ff_ligand.chemcomp_to_rdkit(ml.monomers["TRP"])
        self.assertTrue(any(a.GetIsAromatic() for a in mol.GetAtoms()))
        # subset of atoms (e.g. without OXT/H2/H3 as in a polymer)
        sub = [n for n in names if n not in ("OXT", "H2", "H3")]
        mol2, names2 = ff_ligand.chemcomp_to_rdkit(ml.monomers["TRP"], sub)
        self.assertEqual(mol2.GetNumAtoms(), len(sub))

    def test_ligand_openff_prior_optional(self):
        try:
            import openmm  # noqa: F401
            import openmmforcefields  # noqa: F401
            from openff.toolkit import Molecule  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM/openmmforcefields/openff-toolkit are not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        st = utils.fileio.read_structure("tests/biotin/biotin_talos.pdb")
        monlib = utils.restraints.load_monomer_library(st)
        utils.restraints.add_hydrogens(st, monlib, "nucl")
        self._renumber_serials(st)
        rp = RefineParams(st, refine_xyz=True, adp_mode=0, refine_occ=False, refine_dfrac=False, cfg=None)
        # without ligand parameterisation biotin has no template
        with self.assertRaises(RuntimeError):
            AmberFFPrior(st, rp, platform_name="Reference", monlib=monlib, ligand_ff="none")
        # gasteiger charges keep the test fast and free of AmberTools
        prior = AmberFFPrior(st, rp, platform_name="Reference", monlib=monlib,
                             ligand_ff="openff-2.2.1", ligand_charge="gasteiger")
        self.assertEqual(len(prior.ligand_info), 1)
        self.assertEqual(prior.ligand_info[0][0], "BTN")
        self.assertEqual(prior.ligand_info[0][2], -1)
        ene, grad = prior.calc_target_and_grad()
        self.assertTrue(numpy.isfinite(ene))
        self.assertEqual(len(grad), rp.n_params())
        self.assertTrue(numpy.all(numpy.isfinite(grad)))
        hd = prior.hessian_diag_vector()
        self.assertTrue(numpy.all(hd[hd > 0] >= prior.hessian_diag))
        # gradient check on one atom
        cra = next(c for c in st[0].all() if c.atom.element.name == "C")
        pidx = list(rp.atom_to_param(ext.RefineParams.Type.X))[cra.atom.serial - 1]
        x0 = numpy.array(cra.atom.pos.tolist()); h = 1e-4; num = []
        for ax in range(3):
            es = []
            for sgn in (+1, -1):
                v = x0.copy(); v[ax] += sgn * h; cra.atom.pos = gemmi.Position(*v)
                es.append(prior.calc_target_and_grad(target_only=True)[0])
            num.append((es[0] - es[1]) / (2 * h))
        cra.atom.pos = gemmi.Position(*x0)
        numpy.testing.assert_allclose(num, grad[pidx*3:pidx*3+3], rtol=1e-3, atol=1e-2)


class TestAmberLigandRealData(unittest.TestCase):
    """7db6 (melatonin receptor MT1-Gi1 with ramelteon, EMD-30627, 3.3 A): AMBER prior with an
    OpenFF-parameterised ligand on real cryo-EM half maps. Needs the OpenFF environment and network
    access for the first run (files are cached in tests/7db6, verified by MD5)."""
    root = os.path.abspath(os.path.dirname(__file__))
    url_md5 = (("half1", "https://files.wwpdb.org/pub/emdb/structures/EMD-30627/other/emd_30627_half_map_1.map.gz",
                "319dc338b6a44f9be02f5c4b2c07c1b4"),
               ("half2", "https://files.wwpdb.org/pub/emdb/structures/EMD-30627/other/emd_30627_half_map_2.map.gz",
                "322e28cb5fd70bc893c01c02fbec6a2b"),
               ("mmcif", "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/db/7db6.cif.gz",
                "86fce1c3112c758ff4acc6c0f4bc0e54"))

    @classmethod
    def download(cls):
        import hashlib
        from urllib.request import urlretrieve
        wd = os.path.join(cls.root, "7db6")
        os.makedirs(wd, exist_ok=True)
        data = {}
        for name, url, md5 in cls.url_md5:
            dst = os.path.join(wd, os.path.basename(url))
            if not os.path.exists(dst):
                print("downloading {}".format(url))
                urlretrieve(url, dst)
            with open(dst, "rb") as f:
                md5f = hashlib.md5(f.read()).hexdigest()
            if md5 != md5f:
                raise RuntimeError("md5 mismatch for {}: {} != {}".format(dst, md5f, md5))
            data[name] = dst
        return data

    def test_refine_spa_amber_openff_7db6(self):
        try:
            import openmm  # noqa: F401
            import openmmforcefields  # noqa: F401
            from openff.toolkit import Molecule  # noqa: F401
        except ImportError:
            self.skipTest("OpenMM/openmmforcefields/openff-toolkit are not installed")
        if "CLIBD_MON" not in os.environ:
            self.skipTest("CLIBD_MON is not set")
        import json
        import tempfile
        from servalcat.__main__ import main
        data = self.download()
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory(prefix="servalcat_7db6_") as d:
            os.chdir(d)
            try:
                sys.argv = ["", "refine_spa_norefmac",
                            "--model", data["mmcif"], "--halfmaps", data["half1"], data["half2"],
                            "-d", "3.3", "--ncycle", "2",
                            "--amber_enable", "--amber_his_state", "HIE", "--amber_platform", "Reference",
                            "--hout", "-o", "refined"]
                main()
                self.assertTrue(os.path.isfile("refined.mmcif"))
                self.assertTrue(os.path.isfile("refined_stats.json"))
                with open("refined.log") as ifs:
                    log = ifs.read()
                self.assertIn("JEV parameterised with openff-2.2.1", log)
                self.assertIn("disulfide bond(s) -> CYX", log)
                ff = [float(l.split("=")[1]) for l in log.splitlines() if l.startswith(" ff= ")]
                self.assertGreaterEqual(len(ff), 3)
                self.assertTrue(all(numpy.isfinite(ff)))
                self.assertLess(ff[-1], ff[0])  # AMBER energy must not blow up over the cycles
                stats = json.load(open("refined_stats.json"))
                self.assertEqual(len(stats), 3)  # cycle 0 + 2 cycles
                self.assertTrue(all(s_["fval_decreased"] for s_ in stats[1:]))
                fsc = [float(l.split("=")[1]) for l in log.splitlines() if l.strip().startswith("FSCaverage(full)")]
                self.assertGreater(fsc[-1], 0.6)
                # ligand kept its geometry: JEV heavy atoms and hydrogens still bonded (no atom flew away)
                st = utils.fileio.read_structure("refined.mmcif")
                jev = [r for ch in st[0] for r in ch if r.name == "JEV"]
                self.assertEqual(len(jev), 1)
                self.assertEqual(len(jev[0]), 40)
                ns = gemmi.NeighborSearch(st[0], st.cell, 3).populate()
                for a in jev[0]:
                    dmin = min((m.pos.dist(a.pos) for m in (mk.to_cra(st[0]).atom for mk in ns.find_atoms(a.pos, '\0', radius=2.0))
                                if m is not a), default=9.)
                    self.assertLess(dmin, 1.9)
            finally:
                sys.argv = [""]
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
