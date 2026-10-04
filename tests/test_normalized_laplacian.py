from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.sparse.csgraph import laplacian as scipy_laplacian


REPO_ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"
if str(SYNTHETIC_DIR) not in sys.path:
    sys.path.insert(0, str(SYNTHETIC_DIR))

import _benchmark_utils as benchmark  # noqa: E402

NORMALIZED_DIR = REPO_ROOT / "main_experiment" / "normalized_laplacian"
if str(NORMALIZED_DIR) not in sys.path:
    sys.path.insert(0, str(NORMALIZED_DIR))
import run_normalized_laplacian as experiment  # noqa: E402


class NormalizedOperatorTest(unittest.TestCase):
    def setUp(self) -> None:
        # Weighted, irregular, connected graph with deliberately large self loops.
        # Those loops must not change its Laplacian or normalization degrees.
        self.adjacency = np.array(
            [
                [80.0, 2.0, 0.0, 1.0],
                [2.0, 50.0, 3.0, 0.0],
                [0.0, 3.0, 20.0, 4.0],
                [1.0, 0.0, 4.0, 90.0],
            ]
        )
        self.laplacian = scipy_laplacian(self.adjacency, normed=False)

    def test_normalization_matches_scipy_and_excludes_self_loops(self) -> None:
        actual, degrees = benchmark.normalize_graph_laplacian(self.laplacian)
        expected = scipy_laplacian(self.adjacency, normed=True)
        np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)
        np.testing.assert_array_equal(degrees, [3.0, 5.0, 7.0, 5.0])

        loop_free = self.adjacency.copy()
        np.fill_diagonal(loop_free, 0.0)
        np.testing.assert_allclose(
            actual, scipy_laplacian(loop_free, normed=True), rtol=1e-14, atol=1e-14
        )

    def test_entire_normalized_spectrum_solves_generalized_eigenproblem(self) -> None:
        normalized, degrees = benchmark.normalize_graph_laplacian(self.laplacian)
        eigenvalues, normalized_vectors = eigh(normalized)
        coordinates = normalized_vectors / np.sqrt(degrees)[:, None]
        expected_eigenvalues, _ = eigh(self.laplacian, np.diag(degrees))

        np.testing.assert_allclose(eigenvalues, expected_eigenvalues, atol=1e-14)
        np.testing.assert_allclose(
            self.laplacian @ coordinates,
            degrees[:, None] * coordinates * eigenvalues[None, :],
            atol=1e-14,
        )
        np.testing.assert_allclose(
            coordinates.T @ (degrees[:, None] * coordinates),
            np.eye(len(degrees)),
            atol=1e-14,
        )
        # The normalized null vector is proportional to sqrt(degree),
        # while the generalized/null coordinate is constant.
        np.testing.assert_allclose(normalized @ np.sqrt(degrees), 0.0, atol=1e-14)
        np.testing.assert_allclose(
            coordinates[:, 0], np.full(len(degrees), coordinates[0, 0]), atol=1e-14
        )
        self.assertAlmostEqual(float(degrees @ coordinates[:, 1]), 0.0, places=13)

    def test_invalid_operator_inputs_are_rejected(self) -> None:
        bad_inputs = (
            np.zeros((0, 0)),
            np.ones((2, 3)),
            np.zeros((2, 2)),
            np.array([[1.0, -1.0], [0.0, 1.0]]),
            np.array([[np.nan, -1.0], [-1.0, 1.0]]),
        )
        for operator in bad_inputs:
            with self.subTest(shape=operator.shape):
                with self.assertRaises(ValueError):
                    benchmark.normalize_graph_laplacian(operator)

    def test_regular_two_walk_graph_preserves_combinatorial_fiedler_direction(self) -> None:
        block = np.array([[0.8, 0.2], [0.2, 0.8]])
        laplacian = benchmark.mheatmap_two_walk_laplacian(block, alpha=1.0)
        normalized, degrees = benchmark.normalize_graph_laplacian(laplacian)
        np.testing.assert_allclose(degrees, np.full(4, degrees[0]), atol=1e-14)
        np.testing.assert_allclose(normalized, laplacian / degrees[0], atol=1e-14)
        old_values, old_vectors = eigh(laplacian)
        new_values, new_vectors = eigh(normalized)
        self.assertGreater(float(old_values[2] - old_values[1]), 0.1)
        np.testing.assert_allclose(new_values, old_values / degrees[0], atol=1e-14)
        self.assertAlmostEqual(abs(float(old_vectors[:, 1] @ new_vectors[:, 1])), 1.0)


class NormalizedOrderingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.matrix = np.array(
            [
                [10.0, 2.0, 0.0, 0.0, 0.0],
                [2.0, 9.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        )

    def test_package_calls_retain_original_adjacency_global_scale_and_components(self) -> None:
        original = benchmark.mheatmap_two_walk_laplacian
        for alpha in (0.25, 4.0):
            captured_blocks: list[np.ndarray] = []

            def recording_laplacian(block: np.ndarray, alpha: float) -> np.ndarray:
                captured_blocks.append(np.asarray(block).copy())
                return original(block, alpha=alpha)

            with self.subTest(alpha=alpha):
                with patch.object(
                    benchmark,
                    "mheatmap_two_walk_laplacian",
                    side_effect=recording_laplacian,
                ) as mocked:
                    result = benchmark.normalized_laplacian_reorder(self.matrix, alpha=alpha)

                self.assertEqual(mocked.call_count, 2)
                self.assertEqual([call.kwargs["alpha"] for call in mocked.call_args_list], [alpha, alpha])
                np.testing.assert_array_equal(captured_blocks[0], [[1.0, 0.2], [0.2, 0.9]])
                np.testing.assert_array_equal(captured_blocks[1], [[0.1]])
                self.assertEqual(result.component_count, 2)
                self.assertEqual(result.component_masses, (23.0, 1.0))
                self.assertEqual(len(result.spectral_diagnostics), 2)
                self.assertEqual(set(result.row_order[:2]), {0, 1})
                self.assertEqual(set(result.col_order[:2]), {0, 1})
                np.testing.assert_array_equal(result.row_order[2:], [2, 3])
                np.testing.assert_array_equal(result.col_order[2:], [2, 3, 4])
                for diagnostic in result.spectral_diagnostics:
                    self.assertLess(diagnostic["generalized_relative_residual"], 1e-12)
                    self.assertLess(diagnostic["null_mode_residual"], 1e-12)

    def test_sorting_matches_direct_generalized_eigensolve_for_both_walks(self) -> None:
        # Nonuniform degrees and simple Fiedler eigenvalues avoid ambiguous
        # eigenspace/tie comparisons. scipy.linalg.eigh(L, D) is an independent
        # reference for the declared f = D^-1/2 u ordering convention.
        matrix = np.array([[4.0, 0.0, 1.0], [2.0, 1.0, 0.0], [0.0, 1.0, 3.0]])
        block = matrix / matrix.max()
        for mode in ("one_walk", "two_walk"):
            with self.subTest(mode=mode):
                if mode == "one_walk":
                    adjacency = np.block([[np.zeros((3, 3)), block], [block.T, np.zeros((3, 3))]])
                    laplacian = scipy_laplacian(adjacency, normed=False)
                else:
                    laplacian = benchmark.mheatmap_two_walk_laplacian(block, alpha=4.0)
                degrees = np.diag(laplacian)
                eigenvalues, coordinates = eigh(laplacian, np.diag(degrees))
                self.assertGreater(float(eigenvalues[2] - eigenvalues[1]), 0.01)
                coordinate = coordinates[:, 1]
                if coordinate[int(np.argmax(np.abs(coordinate)))] < 0:
                    coordinate *= -1.0
                joint_order = np.argsort(coordinate, kind="stable")
                expected_rows = joint_order[joint_order < 3]
                expected_cols = joint_order[joint_order >= 3] - 3

                result = benchmark.normalized_laplacian_reorder(matrix, mode=mode, alpha=4.0)
                np.testing.assert_array_equal(result.row_order, expected_rows)
                np.testing.assert_array_equal(result.col_order, expected_cols)
                np.testing.assert_array_equal(result.matrix, matrix[np.ix_(expected_rows, expected_cols)])
                self.assertAlmostEqual(
                    result.spectral_diagnostics[0]["lambda_fiedler"], float(eigenvalues[1]), places=13
                )

    def test_one_walk_uses_no_two_walk_package_call(self) -> None:
        with patch.object(benchmark, "mheatmap_two_walk_laplacian") as mocked:
            result = benchmark.normalized_laplacian_reorder(self.matrix, mode="one_walk")
        mocked.assert_not_called()
        self.assertEqual(result.component_count, 2)
        self.assertEqual(result.component_masses, (23.0, 1.0))
        self.assertEqual(len(result.spectral_diagnostics), 2)

    def test_simple_fiedler_mode_can_still_be_flat_within_both_partitions(self) -> None:
        # Equal complete-bipartite weights make the Fiedler direction separate
        # rows from columns, without defining an ordering within either side.
        # Its simple eigenvalue/gap alone therefore cannot certify seriation.
        result = benchmark.normalized_laplacian_reorder(np.ones((3, 3)), alpha=1.0)
        self.assertEqual(len(result.spectral_diagnostics), 1)
        diagnostic = result.spectral_diagnostics[0]
        self.assertGreater(diagnostic["fiedler_gap"], 0.1)
        self.assertLess(diagnostic["row_relative_spread"], 1e-10)
        self.assertLess(diagnostic["col_relative_spread"], 1e-10)
        self.assertLess(diagnostic["generalized_relative_residual"], 1e-12)

    def test_determinism_and_global_scale_invariance(self) -> None:
        matrix = np.array([[4.0, 0.0, 1.0], [2.0, 1.0, 0.0], [0.0, 1.0, 3.0]])
        for mode in ("one_walk", "two_walk"):
            with self.subTest(mode=mode):
                first = benchmark.normalized_laplacian_reorder(matrix, mode=mode, alpha=0.25)
                repeated = benchmark.normalized_laplacian_reorder(matrix, mode=mode, alpha=0.25)
                scaled = benchmark.normalized_laplacian_reorder(19.0 * matrix, mode=mode, alpha=0.25)
                np.testing.assert_array_equal(first.row_order, repeated.row_order)
                np.testing.assert_array_equal(first.col_order, repeated.col_order)
                np.testing.assert_array_equal(first.row_order, scaled.row_order)
                np.testing.assert_array_equal(first.col_order, scaled.col_order)

    def test_all_zero_matrix_has_identity_order_without_eigensolver(self) -> None:
        matrix = np.zeros((4, 6))
        with patch.object(benchmark, "eigh") as mocked:
            result = benchmark.normalized_laplacian_reorder(matrix)
        mocked.assert_not_called()
        np.testing.assert_array_equal(result.row_order, np.arange(4))
        np.testing.assert_array_equal(result.col_order, np.arange(6))
        np.testing.assert_array_equal(result.matrix, matrix)
        self.assertEqual(result.component_count, 0)
        self.assertEqual(result.spectral_diagnostics, ())

    def test_invalid_alpha_and_mode_are_rejected_including_empty_support(self) -> None:
        for matrix in (self.matrix, np.zeros((2, 3))):
            for alpha in (0.0, -1.0, float("nan"), float("inf")):
                with self.subTest(alpha=alpha, matrix_sum=float(matrix.sum())):
                    with self.assertRaisesRegex(ValueError, "positive and finite"):
                        benchmark.normalized_laplacian_reorder(matrix, alpha=alpha)
            with self.assertRaisesRegex(ValueError, "unknown spectral mode"):
                benchmark.normalized_laplacian_reorder(matrix, mode="unknown")


class NormalizedAnalysisTest(unittest.TestCase):
    def test_adaptive_selection_uses_table1_tolerance_and_smaller_alpha_policy(self) -> None:
        # An exact maximum would pick alpha=2 in the first group, but the
        # submitted Table 1 tolerance deliberately retains alpha=1.
        for increment, expected_alpha in ((0.0, 1.0), (0.5e-12, 1.0), (1.5e-12, 2.0)):
            with self.subTest(increment=increment):
                group = pd.DataFrame({"alpha": [2.0, 1.0], "mwb_auc": [0.8 + increment, 0.8]})
                self.assertEqual(float(experiment.select_adaptive(group)["alpha"]), expected_alpha)
        # The policy compares each candidate to the retained score, not to
        # a separately computed maximum: near ties cannot silently accumulate.
        sequential = pd.DataFrame({
            "alpha": [4.0, 2.0, 1.0], "mwb_auc": [0.8 + 1.8e-12, 0.8 + 0.9e-12, 0.8]
        })
        self.assertEqual(float(experiment.select_adaptive(sequential)["alpha"]), 4.0)
        with self.assertRaisesRegex(ValueError, "empty group"):
            experiment.select_adaptive(pd.DataFrame(columns=["alpha", "mwb_auc"]))

    def test_protocols_keep_fixed_alphas_and_select_each_variant_independently(self) -> None:
        metadata = {
            "suite": "synthetic", "case_id": "case-a", "family_key": "family-a",
            "family": "Family A", "size_key": "small", "size": "Small", "dataset": "",
            "seed": 0, "n_rows": 3, "n_cols": 3, "matrix_sha256": "fixture",
        }
        rows = []
        for laplacian, best_alpha, offset in (("combinatorial", 2.0, 0.0), ("normalized", 8.0, 0.1)):
            rows.append({
                **metadata, "laplacian": laplacian, "method": "OW", "alpha": np.nan,
                "two_sum": 0.5 + offset, "band_mass_10": 0.6, "mwb_auc": 0.7,
            })
            for alpha in experiment.ALPHAS:
                rows.append({
                    **metadata, "laplacian": laplacian, "method": "TW", "alpha": alpha,
                    "two_sum": alpha / 100.0 + offset, "band_mass_10": alpha / 20.0,
                    "mwb_auc": 0.9 if alpha == best_alpha else 0.4,
                })
        protocols = experiment.protocol_rows(pd.DataFrame(rows)).set_index("variant")
        self.assertEqual(set(protocols.index), set(experiment.VARIANTS))
        for prefix in ("c", "n"):
            self.assertTrue(np.isnan(protocols.loc[f"{prefix}_ow", "selected_alpha"]))
            self.assertEqual(float(protocols.loc[f"{prefix}_tw_fixed_1", "selected_alpha"]), 1.0)
            self.assertEqual(float(protocols.loc[f"{prefix}_tw_fixed_12", "selected_alpha"]), 12.0)
        self.assertEqual(float(protocols.loc["c_tw_adaptive", "selected_alpha"]), 2.0)
        self.assertEqual(float(protocols.loc["n_tw_adaptive", "selected_alpha"]), 8.0)
        self.assertAlmostEqual(float(protocols.loc["c_tw_adaptive", "two_sum"]), 0.02)
        self.assertAlmostEqual(float(protocols.loc["n_tw_adaptive", "two_sum"]), 0.18)

    @staticmethod
    def _paired_fixture() -> pd.DataFrame:
        rows = []
        groups = (
            ("synthetic", "Family A", "Small", "", (0.125, 0.25, 0.375), 2),
            ("synthetic", "Family B", "Large", "", (0.375, 0.5, 0.125), 2),
            ("real", "", "", "Dataset A", (0.25, -0.25, 0.125), 1),
        )
        for suite, family, size, dataset, differences, count in groups:
            for seed in range(count):
                for _, control, normalized in experiment.COMPARISONS:
                    for variant in (control, normalized):
                        values = np.full(3, 0.5)
                        if variant == normalized:
                            values += differences
                        rows.append({
                            "suite": suite, "case_id": f"{suite}/{family}/{size}/{dataset}/{seed}",
                            "family": family, "size": size, "dataset": dataset,
                            "seed": seed,
                            "variant": variant,
                            **dict(zip(experiment.METRICS, values, strict=True)),
                        })
        return pd.DataFrame(rows)

    def test_paired_bootstrap_constant_regimes_have_exact_intervals_and_stratified_mean(self) -> None:
        # Every seed within a regime has the same paired delta. Therefore
        # each bootstrap draw must reproduce the group delta exactly, and
        # the overall stratified draw is the equally weighted regime mean.
        # Equal seed counts mirror the actual 20-seed/regime experiment.
        results = experiment.paired_deltas(self._paired_fixture())
        synthetic = results.loc[results["suite"] == "synthetic"]
        groups = synthetic.loc[synthetic["group_scope"] == "group"]
        expected = {
            "Family A": dict(zip(experiment.METRICS, (0.125, 0.25, 0.375), strict=True)),
            "Family B": dict(zip(experiment.METRICS, (0.375, 0.5, 0.125), strict=True)),
        }
        for row in groups.itertuples():
            delta = expected[row.family][row.metric]
            self.assertEqual(row.delta_mean, delta)
            self.assertEqual(row.delta_std, 0.0)
            self.assertEqual(row.delta_ci_low, delta)
            self.assertEqual(row.delta_ci_high, delta)
            self.assertEqual(row.n_cases, 2)
        overall = synthetic.loc[synthetic["group_scope"] == "overall"]
        expected_overall = dict(zip(experiment.METRICS, (0.25, 0.375, 0.25), strict=True))
        for row in overall.itertuples():
            self.assertEqual(row.delta_mean, expected_overall[row.metric])
            self.assertEqual(row.delta_ci_low, expected_overall[row.metric])
            self.assertEqual(row.delta_ci_high, expected_overall[row.metric])
            self.assertEqual(row.n_cases, 4)
            self.assertEqual(row.normalized_win_rate, 0.0 if row.metric == "two_sum" else 1.0)
            self.assertEqual(row.tie_rate, 0.0)

    def test_overall_bootstrap_resamples_matching_seeds_jointly_across_regimes(self) -> None:
        # Regimes use the same generator seeds. Opposite seed-matched effects
        # cancel in every macro-average. Independently resampling regimes
        # would invent variability and produce a nonzero overall CI width.
        rows = []
        seed_effects = ((0.125, 0.25, 0.125), (0.375, 0.125, 0.25))
        for family, sign in (("Family A", 1.0), ("Family B", -1.0)):
            for seed, effects in enumerate(seed_effects):
                for _, control, normalized in experiment.COMPARISONS:
                    for variant in (control, normalized):
                        values = np.full(3, 0.5)
                        if variant == normalized:
                            values += sign * np.asarray(effects)
                        rows.append({
                            "suite": "synthetic", "case_id": f"{family}/seed-{seed}",
                            "family": family, "size": "Small", "dataset": "", "seed": seed,
                            "variant": variant,
                            **dict(zip(experiment.METRICS, values, strict=True)),
                        })
        results = experiment.paired_deltas(pd.DataFrame(rows))
        groups = results.loc[results["group_scope"] == "group"]
        self.assertTrue((groups["delta_std"] > 0.0).all())
        self.assertTrue((groups["delta_ci_high"] > groups["delta_ci_low"]).all())
        overall = results.loc[results["group_scope"] == "overall"]
        self.assertEqual(len(overall), len(experiment.COMPARISONS) * len(experiment.METRICS))
        np.testing.assert_array_equal(overall["delta_mean"], 0.0)
        np.testing.assert_array_equal(overall["delta_ci_low"], 0.0)
        np.testing.assert_array_equal(overall["delta_ci_high"], 0.0)

    def test_real_datasets_report_paired_deltas_without_seed_confidence_intervals(self) -> None:
        results = experiment.paired_deltas(self._paired_fixture())
        real = results.loc[results["suite"] == "real"]
        self.assertEqual(len(real), len(experiment.COMPARISONS) * len(experiment.METRICS) * 2)
        self.assertTrue(real[["delta_ci_low", "delta_ci_high", "delta_std"]].isna().all().all())
        expected = dict(zip(experiment.METRICS, (0.25, -0.25, 0.125), strict=True))
        for row in real.itertuples():
            self.assertEqual(row.delta_mean, expected[row.metric])
            self.assertEqual(row.n_cases, 1)
            self.assertEqual(row.normalized_win_rate, 1.0 if row.metric == "mwb_auc" else 0.0)


if __name__ == "__main__":
    unittest.main()
