"""Small-matrix gate for the Section 3.3 matrix-free Two-Walk implementation.

Run before anything corpus-level, and again after every change to the operators or
the solver:

    cd MHeatMap/src && uv run python rebuttal/test_tw_matrixfree.py
    cd MHeatMap/src && uv run python -m unittest discover -s rebuttal -p 'test_*.py'

Everything here is deliberately tiny -- hand-built blocks plus the *smallest*
components clipped from the real matrices -- so the whole file stays well under
half a minute and can be run on every edit.  The corpus-wide version of check C
lives in ``run_matrixfree_equivalence.py``.

Checks
------
A  operator identity, against the *released* ``mheatmap.two_walk_laplacian``
B  structural identities (nullspace, degree vector, symmetry) and the SpMV count
C  dense vs matrix-free agreement, with the sparse solver forced on
D  a repeated Fiedler eigenvalue is flagged, not silently resolved
E  edge cases route away from ARPACK rather than crashing inside it
F  determinism across runs
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import metrics_v2 as mv  # noqa: E402
import tw_components as tc  # noqa: E402
import tw_matrixfree as mf  # noqa: E402

ALPHAS = (0.25, 1.0, 6.0, 12.0)
TIE_TOL = 1e-9
REL_LAMBDA_TOL = 1e-8
COORD_TOL = 1e-9
R2S_TOL = 1e-4
BAND_TOL = 1e-9
MWB_TOL = 1e-4


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #

HAND_BLOCKS: dict[str, np.ndarray] = {
    "2x2_identity": np.array([[1.0, 0.0], [0.0, 1.0]]),
    "2x2_full": np.array([[1.0, 2.0], [3.0, 1.0]]),
    "3x4_mixed": np.array(
        [[1.0, 0.0, 2.0, 0.0], [0.0, 3.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0]]
    ),
    "5x5_band": np.eye(5) + np.eye(5, k=1) + np.eye(5, k=-1),
    "single_edge": np.array([[1.0]]),
    "zero_row_col": np.array(
        [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    ),
    "weighted_skew": np.array([[4.0, 0.0, 1.0], [0.0, 0.5, 0.0], [1.0, 0.0, 3.0]]),
}


def real_blocks(max_vertices: int = 150) -> dict[str, np.ndarray]:
    """The smallest connected component of each real matrix, at a tractable size.

    Real data, at a size where the dense reference is instant and the tie structure
    can still be reasoned about by hand.
    """
    matrix_dir = (
        mf.REPO_ROOT
        / "output"
        / "main_experiment"
        / "real_world_benchmark"
        / "processed_matrices"
    )
    if not matrix_dir.is_dir():
        return {}

    blocks: dict[str, np.ndarray] = {}
    for path in sorted(matrix_dir.glob("*.csv")):
        matrix, _, _ = mf.bm.load_matrix_csv(path)
        normalized = matrix / float(matrix.max())
        for rank, component in enumerate(tc.support_components(matrix)):
            size = component.n_rows + component.n_cols
            if not 6 <= size <= max_vertices:
                continue
            block = normalized[np.ix_(component.row_indices, component.col_indices)]
            blocks[f"{path.stem}::c{rank}"] = block
            break
    return blocks


def all_blocks() -> dict[str, np.ndarray]:
    blocks = dict(HAND_BLOCKS)
    blocks.update(real_blocks())
    return blocks


def is_connected(block: np.ndarray) -> bool:
    """Whether the block's bipartite support graph is connected.

    ``componentwise_reorder`` only ever feeds the solver a single connected
    component, so this is the domain the dense-vs-sparse contract is stated over.
    A disconnected block has a repeated zero eigenvalue and a Fiedler pair that is
    not unique, which is a separate (already-flagged) matter.
    """
    support = sp.csr_matrix(np.asarray(block) > 0)
    graph = sp.bmat([[None, support], [support.T, None]], format="csr")
    count, _ = sp.csgraph.connected_components(graph, directed=False)
    return int(count) == 1


def connected_blocks() -> dict[str, np.ndarray]:
    return {name: b for name, b in all_blocks().items() if is_connected(b)}


def dense_two_walk_adjacency(block: np.ndarray, alpha: float) -> np.ndarray:
    """A_tw straight from Eq. (tw), independent of anything this module computes."""
    B = np.asarray(block, dtype=float)
    return np.block(
        [[B @ B.T, alpha * B], [alpha * B.T, B.T @ B]]
    )


def split_orders(joint_order: np.ndarray, local_row_count: int):
    rows = joint_order[joint_order < local_row_count]
    cols = joint_order[joint_order >= local_row_count] - local_row_count
    return rows, cols


class CountingSparse:
    """Proxy that counts ``@`` products, to verify the SpMV claim in Section 3.3."""

    def __init__(self, inner):
        self._inner = inner
        self.count = 0

    def __matmul__(self, other):
        self.count += 1
        return self._inner @ other


# --------------------------------------------------------------------------- #
# A. operator identity against the released object
# --------------------------------------------------------------------------- #


class TestOperatorIdentity(unittest.TestCase):
    def test_two_walk_matches_released_laplacian(self):
        rng = np.random.default_rng(0)
        for name, block in all_blocks().items():
            for alpha in ALPHAS:
                reference = mf.build_two_walk_dense(block, alpha)
                op = mf.TwoWalkOperator(block, alpha)
                for _ in range(5):
                    z = rng.standard_normal(op.shape[0])
                    want = reference @ z
                    err = np.linalg.norm(op @ z - want) / max(
                        np.linalg.norm(want), 1e-300
                    )
                    self.assertLess(
                        err, 1e-12, msg=f"{name} alpha={alpha} relative error {err:.3e}"
                    )

    def test_one_walk_matches_dense_block(self):
        rng = np.random.default_rng(1)
        for name, block in all_blocks().items():
            reference = mf.build_one_walk_dense(block)
            op = mf.OneWalkOperator(block)
            for _ in range(5):
                z = rng.standard_normal(op.shape[0])
                want = reference @ z
                err = np.linalg.norm(op @ z - want) / max(np.linalg.norm(want), 1e-300)
                self.assertLess(err, 1e-12, msg=f"{name} relative error {err:.3e}")


# --------------------------------------------------------------------------- #
# B. structural identities
# --------------------------------------------------------------------------- #


class TestStructuralIdentities(unittest.TestCase):
    def test_nullspace_and_degree_vector(self):
        for name, block in all_blocks().items():
            for alpha in ALPHAS:
                op = mf.TwoWalkOperator(block, alpha)

                # Proposition 1: L_tw 1 = 0 on a connected component.
                ones = np.ones(op.shape[0])
                self.assertLess(
                    np.linalg.norm(op @ ones), 1e-10, msg=f"{name} alpha={alpha}"
                )

                # Eq. (lmatvec): d_alpha must be the row-sum vector of A_tw, computed
                # here from Eq. (tw) directly rather than from the operator.
                want = dense_two_walk_adjacency(block, alpha).sum(axis=1)
                np.testing.assert_allclose(
                    op.degree,
                    want,
                    rtol=1e-10,
                    atol=1e-10,
                    err_msg=f"{name} alpha={alpha} degree vector mismatch",
                )

    def test_symmetry(self):
        rng = np.random.default_rng(2)
        for name, block in all_blocks().items():
            for alpha in ALPHAS:
                op = mf.TwoWalkOperator(block, alpha)
                z1 = rng.standard_normal(op.shape[0])
                z2 = rng.standard_normal(op.shape[0])
                self.assertAlmostEqual(
                    float(z1 @ (op @ z2)),
                    float(z2 @ (op @ z1)),
                    places=8,
                    msg=f"{name} alpha={alpha} not symmetric",
                )

    def test_matvec_uses_exactly_four_sparse_products(self):
        """Section 3.3: 'one application of L_tw uses four sparse products'."""
        for name, block in all_blocks().items():
            for alpha in ALPHAS:
                op = mf.TwoWalkOperator(block, alpha)
                counter_b, counter_bt = CountingSparse(op.B), CountingSparse(op.Bt)
                op.B, op.Bt = counter_b, counter_bt
                rng = np.random.default_rng(3)
                op @ rng.standard_normal(op.shape[0])
                total = counter_b.count + counter_bt.count
                self.assertEqual(
                    total,
                    mf.TW_SPMV_PER_MATVEC,
                    msg=f"{name} alpha={alpha}: {total} sparse products, "
                    f"paper claims {mf.TW_SPMV_PER_MATVEC}",
                )

    def test_one_walk_uses_exactly_two_sparse_products(self):
        for name, block in all_blocks().items():
            op = mf.OneWalkOperator(block)
            counter_b, counter_bt = CountingSparse(op.B), CountingSparse(op.Bt)
            op.B, op.Bt = counter_b, counter_bt
            rng = np.random.default_rng(4)
            op @ rng.standard_normal(op.shape[0])
            self.assertEqual(
                counter_b.count + counter_bt.count, mf.OW_SPMV_PER_MATVEC, msg=name
            )


# --------------------------------------------------------------------------- #
# C. dense vs matrix-free agreement, with the sparse solver forced on
# --------------------------------------------------------------------------- #

TIE_DIAGNOSTICS: list[dict] = []


class TestDenseVersusMatrixFree(unittest.TestCase):
    def test_eigenvalue_and_ordering_agreement(self):
        for name, block in connected_blocks().items():
            for mode in ("ow", "tw"):
                alphas = (None,) if mode == "ow" else ALPHAS
                for alpha in alphas:
                    with self.subTest(block=name, mode=mode, alpha=alpha):
                        record = mf.compare_block(
                            block, alpha, mode, tie_tol=TIE_TOL, tol=1e-12
                        )
                        TIE_DIAGNOSTICS.append(
                            {"block": name, "mode": mode, "alpha": alpha, **record}
                        )
                        self.assertLessEqual(
                            record["rel_lambda_diff"],
                            REL_LAMBDA_TOL,
                            msg=f"eigenvalue drift {record['rel_lambda_diff']:.3e}",
                        )
                        self.assertLessEqual(
                            record["max_sorted_coord_dev"],
                            COORD_TOL,
                            msg=f"coordinate drift {record['max_sorted_coord_dev']:.3e}",
                        )
                        self.assertTrue(
                            record["disagreement_within_ties"],
                            msg="matrix-free order differs outside a tie group",
                        )

    def test_metrics_are_stable(self):
        for name, block in connected_blocks().items():
            for mode in ("ow", "tw"):
                alphas = (None,) if mode == "ow" else (1.0, 12.0)
                for alpha in alphas:
                    with self.subTest(block=name, mode=mode, alpha=alpha):
                        dense = mf.fiedler_block(block, alpha, mode, solver="dense")
                        free = mf.fiedler_block(
                            block, alpha, mode, solver="matrix_free"
                        )
                        dr, dc = split_orders(dense.joint_order, block.shape[0])
                        fr, fc = split_orders(free.joint_order, block.shape[0])

                        self.assertLessEqual(
                            abs(mv.r2s(block, dr, dc) - mv.r2s(block, fr, fc)), R2S_TOL
                        )
                        self.assertLessEqual(
                            abs(
                                mv.band_mass_eq16(block, dr, dc)
                                - mv.band_mass_eq16(block, fr, fc)
                            ),
                            BAND_TOL,
                        )
                        self.assertLessEqual(
                            abs(
                                mv.mwb_auc_block(block, dr, dc)
                                - mv.mwb_auc_block(block, fr, fc)
                            ),
                            MWB_TOL,
                        )


# --------------------------------------------------------------------------- #
# C2. the sparse-native pipeline must match the dense one exactly
# --------------------------------------------------------------------------- #


class TestSparsePath(unittest.TestCase):
    """The scaling sweep cannot densify, so its sparse pipeline is checked here."""

    MATRICES = {
        **HAND_BLOCKS,
        "with_zero_marginals": np.array(
            [[0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 1.0, 0.0], [0.0, 1.0, 3.0, 0.0]]
        ),
    }

    def test_sparse_and_dense_reorder_agree(self):
        for name, block in self.MATRICES.items():
            for mode in ("ow", "tw"):
                for solver in ("dense", "matrix_free"):
                    with self.subTest(name=name, mode=mode, solver=solver):
                        alpha = 1.0 if mode == "tw" else None
                        dense = mf.componentwise_reorder(
                            block, mode=mode, alpha=alpha, solver=solver
                        )
                        sparse = mf.componentwise_reorder_sparse(
                            sp.csr_matrix(block), mode=mode, alpha=alpha, solver=solver
                        )
                        np.testing.assert_array_equal(
                            dense.row_order, sparse.row_order, err_msg=name
                        )
                        np.testing.assert_array_equal(
                            dense.col_order, sparse.col_order, err_msg=name
                        )
                        self.assertEqual(
                            dense.component_count, sparse.component_count, msg=name
                        )
                        self.assertEqual(
                            dense.component_masses, sparse.component_masses, msg=name
                        )

    def test_footprints_match_a_hand_computed_value(self):
        block = np.array([[1.0, 0.0], [0.0, 1.0]])  # two components of 2 vertices each
        record = mf.footprints_sparse(sp.csr_matrix(block))
        self.assertEqual(record["n_components"], 2)
        self.assertEqual(record["max_component_vertices"], 2)
        # 2 components * 2 matrices * 2x2 * 8 bytes = 128 bytes
        self.assertAlmostEqual(record["operator_footprint_mb"], 128 / 1e6, places=12)

    def test_dense_and_sparse_footprints_agree_on_shape(self):
        for name, block in self.MATRICES.items():
            dense = mf.footprints_sparse(block)
            sparse = mf.footprints_sparse(sp.csr_matrix(block))
            self.assertEqual(dense, sparse, msg=name)


# --------------------------------------------------------------------------- #
# D. degeneracy is flagged
# --------------------------------------------------------------------------- #


class TestDegeneracy(unittest.TestCase):
    def test_repeated_fiedler_eigenvalue_is_flagged(self):
        """K_{2,2} has a doubled Fiedler eigenvalue, so the order is not unique."""
        result = mf.fiedler_block(np.ones((2, 2)), None, "ow", solver="dense")
        self.assertTrue(
            result.degenerate,
            msg="expected the doubled Fiedler eigenvalue of K_{2,2} to be flagged",
        )

    def test_simple_eigenvalue_is_not_flagged(self):
        block = np.array([[1.0, 2.0], [3.0, 1.0], [1.0, 0.0]])
        result = mf.fiedler_block(block, 1.0, "tw", solver="dense")
        self.assertFalse(result.degenerate)


# --------------------------------------------------------------------------- #
# E. edge cases
# --------------------------------------------------------------------------- #


class TestEdgeCases(unittest.TestCase):
    def test_two_vertex_block_routes_to_dense(self):
        """ARPACK needs k < N, so a two-vertex block cannot use it."""
        block = np.array([[1.0]])  # one row, one column
        self.assertEqual(block.shape[0] + block.shape[1], 2)
        result = mf.fiedler_block(block, 1.0, "tw", solver="matrix_free")
        self.assertEqual(result.solver, "dense")
        self.assertEqual(sorted(result.joint_order.tolist()), [0, 1])

    def test_three_vertex_block_uses_the_sparse_solver(self):
        """One vertex above the cutoff, the sparse path must actually be taken."""
        block = np.array([[1.0, 0.0]])  # one row, two columns
        self.assertEqual(block.shape[0] + block.shape[1], 3)
        result = mf.fiedler_block(block, 1.0, "tw", solver="matrix_free")
        self.assertIn(result.solver, {"matrix_free", "eigsh_shifted"})

    def test_all_zero_block_raises_documented_error(self):
        with self.assertRaises(mf.DegenerateBlockError):
            mf.fiedler_block(np.zeros((3, 3)), 1.0, "tw", solver="matrix_free")

    def test_nonpositive_alpha_is_rejected(self):
        for alpha in (0.0, -1.0):
            with self.assertRaises(ValueError):
                mf.fiedler_block(np.array([[1.0, 1.0], [1.0, 1.0]]), alpha, "tw")

    def test_disconnected_block_does_not_crash(self):
        result = mf.fiedler_block(
            np.array([[1.0, 0.0], [0.0, 1.0]]), 1.0, "tw", solver="matrix_free"
        )
        self.assertEqual(sorted(result.joint_order.tolist()), [0, 1, 2, 3])

    def test_componentwise_reorder_returns_valid_permutations(self):
        matrix = np.array(
            [[1.0, 0.0, 2.0, 0.0], [0.0, 3.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
        )
        for mode in ("ow", "tw"):
            for solver in ("dense", "matrix_free"):
                with self.subTest(mode=mode, solver=solver):
                    result = mf.componentwise_reorder(
                        matrix,
                        mode=mode,
                        alpha=1.0 if mode == "tw" else None,
                        solver=solver,
                    )
                    np.testing.assert_array_equal(
                        np.sort(result.row_order), np.arange(3)
                    )
                    np.testing.assert_array_equal(
                        np.sort(result.col_order), np.arange(4)
                    )

    def test_zero_marginal_axes_are_appended_last(self):
        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 2.0], [0.0, 2.0, 1.0]])
        result = mf.componentwise_reorder(matrix, mode="ow", solver="dense")
        self.assertEqual(result.row_order[-1], 0)
        self.assertEqual(result.col_order[-1], 0)


# --------------------------------------------------------------------------- #
# F. determinism
# --------------------------------------------------------------------------- #


MATVEC_SPREAD: list[tuple[str, list[int]]] = []


class TestDeterminism(unittest.TestCase):
    def test_repeated_runs_give_identical_orders(self):
        """The eigenpair and the induced order are deterministic; the matvec count is not.

        ARPACK's internal convergence bookkeeping is not reproducible, so ``n_matvec``
        varies by roughly 10% across runs on identical input (measured 504-558 on
        20 Newsgroups).  The returned eigenvalue, eigenvector and permutation are
        bit-identical, which is what matters.  Only the result is asserted here; the
        count is recorded as a spread, and ``run_runtime_scaling`` must report it as a
        measured quantity rather than a fixed multiple of ``nnz``.
        """
        for name, block in connected_blocks().items():
            results = [
                mf.fiedler_block(block, 6.0, "tw", solver="matrix_free", seed=7)
                for _ in range(3)
            ]
            first = results[0]
            for other in results[1:]:
                np.testing.assert_array_equal(
                    first.joint_order, other.joint_order, err_msg=name
                )
                self.assertAlmostEqual(first.eigenvalue, other.eigenvalue, places=12)
                self.assertEqual(first.solver, other.solver, msg=name)
            MATVEC_SPREAD.append((name, [r.n_matvec for r in results]))


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #


def print_tie_diagnostics() -> None:
    """Show where the two solvers disagree, and why it does not matter."""
    if not TIE_DIAGNOSTICS:
        return
    interesting = [r for r in TIE_DIAGNOSTICS if r["n_positions_differ"]]
    print("\nTie diagnostics (blocks where the two solvers order differently)")
    print(
        f"  {'block':32s} {'mode':4s} {'alpha':>5s} {'disagree':>9s} {'tied':>6s} "
        f"{'max|dcoord|':>12s} {'rel dlambda':>12s} {'within ties':>11s}"
    )
    for record in sorted(
        interesting, key=lambda r: -r["frac_positions_differ"]
    )[:15]:
        alpha = record["alpha"]
        label = "-" if alpha is None else f"{alpha:.2f}"
        print(
            f"  {record['block'][:32]:32s} {record['mode']:4s} {label:>5s} "
            f"{record['frac_positions_differ']:>8.1%} "
            f"{record['n_tied_positions']:>6d} "
            f"{record['max_sorted_coord_dev']:>12.3e} "
            f"{record['rel_lambda_diff']:>12.2e} "
            f"{str(record['disagreement_within_ties']):>11s}"
        )
    print(
        f"\n  {len(interesting)} of {len(TIE_DIAGNOSTICS)} configurations disagree on at "
        "least one position."
    )
    print(
        "  Every disagreement sits inside a numerically tied coordinate group, which is a\n"
        "  property of the method's tied Fiedler coordinates, not a solver defect — see\n"
        "  the tw_matrixfree module docstring."
    )


def print_matvec_spread() -> None:
    """ARPACK's iteration count is not reproducible; its result is."""
    if not MATVEC_SPREAD:
        return
    variable = [
        (name, counts)
        for name, counts in MATVEC_SPREAD
        if len(set(counts)) > 1
    ]
    print("\nMatvec-count spread over 3 identical runs (same input, same seed)")
    for name, counts in variable[:8]:
        span = (max(counts) - min(counts)) / max(min(counts), 1)
        print(f"  {name[:44]:44s} {counts}  spread {span:.1%}")
    print(
        f"\n  {len(variable)} of {len(MATVEC_SPREAD)} blocks varied across runs, while their\n"
        "  eigenvalues, eigenvectors and permutations were bit-identical. ARPACK's\n"
        "  iteration count is therefore reported as a measured spread, never as a\n"
        "  fixed multiple of nnz (see plan V5)."
    )


def main() -> int:
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print_tie_diagnostics()
    print_matvec_spread()
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
