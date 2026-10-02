from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"
REAL_DIR = REPO_ROOT / "main_experiment" / "real_world_benchmark"
for path in (SYNTHETIC_DIR, REAL_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import _benchmark_utils as benchmark  # noqa: E402
from _benchmark_utils import (  # noqa: E402
    DEFAULT_WIDTH_GRID,
    legacy_integer_band_mass,
    load_matrix_csv,
    mwb_auc,
    normalized_band_mass,
    one_walk_reorder,
    support_component_count,
    tw_alpha_reorder,
    tw_reorder_for_alpha,
)
from run_synthetic_evaluation import FAMILIES, SIZES, build_case  # noqa: E402

ALPHA_DIR = REPO_ROOT / "main_experiment" / "alpha_sensitivity"
if str(ALPHA_DIR) not in sys.path:
    sys.path.insert(0, str(ALPHA_DIR))
from run_alpha_sensitivity import (  # noqa: E402
    alpha_zero_diagnostic,
    choose_global_alpha,
    select_best_rows,
)


MATRIX_DIR = (
    REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "processed_matrices"
)


class MetricDefinitionTest(unittest.TestCase):
    def test_equation_16_matches_direct_normalized_coordinate_mask(self) -> None:
        matrix = np.arange(1.0, 1.0 + 8 * 125).reshape(8, 125)
        rows = np.arange(8, dtype=float)[:, None] / 7.0
        cols = np.arange(125, dtype=float)[None, :] / 124.0
        expected = float(matrix[np.abs(rows - cols) <= 0.10].sum() / matrix.sum())
        self.assertAlmostEqual(normalized_band_mass(matrix, 0.10), expected, places=15)

    def test_mwb_is_exact_legacy_band_integral(self) -> None:
        rng = np.random.default_rng(7)
        matrix = rng.integers(0, 8, size=(17, 9)).astype(float)
        scores = np.array(
            [legacy_integer_band_mass(matrix, float(w)) for w in DEFAULT_WIDTH_GRID]
        )
        expected = float(
            np.trapezoid(scores, DEFAULT_WIDTH_GRID)
            / (DEFAULT_WIDTH_GRID[-1] - DEFAULT_WIDTH_GRID[0])
        )
        self.assertEqual(mwb_auc(matrix, DEFAULT_WIDTH_GRID), expected)

    def test_mbta_anchor_reconciles_submitted_and_equation_16_scores(self) -> None:
        matrix, _, _ = load_matrix_csv(MATRIX_DIR / "gtfs_clean_original_matrix.csv")
        self.assertAlmostEqual(legacy_integer_band_mass(matrix, 0.10), 0.023745238, places=8)
        self.assertAlmostEqual(normalized_band_mass(matrix, 0.10), 0.227278504, places=8)


class ComponentwiseOrderingTest(unittest.TestCase):
    def setUp(self) -> None:
        # Two active components with different scales, plus a zero row/column.
        self.matrix = np.array(
            [
                [10.0, 2.0, 0.0, 0.0, 0.0],
                [2.0, 10.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        )

    def test_tw_calls_mheatmap_once_per_component_for_every_alpha(self) -> None:
        original = benchmark.mheatmap_two_walk_laplacian
        for alpha in (1.0, 4.0):
            captured_blocks: list[np.ndarray] = []

            def recording_laplacian(block: np.ndarray, alpha: float) -> np.ndarray:
                captured_blocks.append(np.asarray(block).copy())
                return original(block, alpha=alpha)

            with patch.object(
                benchmark,
                "mheatmap_two_walk_laplacian",
                side_effect=recording_laplacian,
            ) as mocked:
                result = tw_reorder_for_alpha(self.matrix, alpha=alpha)

            self.assertEqual(mocked.call_count, 2)
            self.assertEqual(result.component_count, 2)
            self.assertEqual(result.component_masses, (24.0, 1.0))
            self.assertEqual([float(block.max()) for block in captured_blocks], [1.0, 0.1])
            self.assertEqual(set(result.row_order[:2]), {0, 1})
            self.assertEqual(int(result.row_order[2]), 2)
            self.assertEqual(set(result.col_order[:2]), {0, 1})
            self.assertEqual(int(result.col_order[2]), 2)
            self.assertEqual(int(result.row_order[-1]), 3)
            self.assertEqual(int(result.col_order[-1]), 4)

    def test_one_walk_uses_same_component_and_zero_axis_policy(self) -> None:
        with patch.object(benchmark, "mheatmap_two_walk_laplacian") as mocked:
            result = one_walk_reorder(self.matrix)
        mocked.assert_not_called()
        self.assertEqual(result.component_count, 2)
        self.assertEqual(result.component_masses, (24.0, 1.0))
        self.assertEqual(set(result.row_order[:2]), {0, 1})
        self.assertEqual(int(result.row_order[2]), 2)
        self.assertEqual(set(result.col_order[:2]), {0, 1})
        self.assertEqual(int(result.col_order[2]), 2)
        self.assertEqual(int(result.row_order[-1]), 3)
        self.assertEqual(int(result.col_order[-1]), 4)

    def test_scale_invariance_and_determinism(self) -> None:
        for reorder in (
            one_walk_reorder,
            lambda matrix: tw_alpha_reorder(matrix, alpha=6.0),
        ):
            first = reorder(self.matrix)
            repeated = reorder(self.matrix)
            scaled = reorder(19.0 * self.matrix)
            np.testing.assert_array_equal(first.row_order, repeated.row_order)
            np.testing.assert_array_equal(first.col_order, repeated.col_order)
            np.testing.assert_array_equal(first.row_order, scaled.row_order)
            np.testing.assert_array_equal(first.col_order, scaled.col_order)

    def test_mheatmap_laplacian_is_exact_equation_5(self) -> None:
        block = np.array([[0.2, 1.0], [0.6, 0.4], [0.0, 0.8]])
        alpha = 6.0
        adjacency = np.block(
            [
                [block @ block.T, alpha * block],
                [alpha * block.T, block.T @ block],
            ]
        )
        expected = np.diag(adjacency.sum(axis=1)) - adjacency
        actual = benchmark.mheatmap_two_walk_laplacian(block, alpha=alpha)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-15)

    def test_positive_subunit_alpha_is_supported_and_deterministic(self) -> None:
        first = tw_alpha_reorder(self.matrix, alpha=0.25)
        repeated = tw_alpha_reorder(self.matrix, alpha=0.25)
        np.testing.assert_array_equal(first.row_order, repeated.row_order)
        np.testing.assert_array_equal(first.col_order, repeated.col_order)

    def test_alpha_zero_is_rejected_as_an_order_but_has_two_modes_per_component(self) -> None:
        with self.assertRaisesRegex(ValueError, "alpha must be positive"):
            tw_alpha_reorder(self.matrix, alpha=0.0)
        diagnostic = alpha_zero_diagnostic(self.matrix, {"case_id": "unit"})
        self.assertEqual(diagnostic["component_count"], 2)
        self.assertEqual(diagnostic["expected_nullity"], 4)
        self.assertEqual(diagnostic["observed_nullity"], 4)
        self.assertFalse(diagnostic["joint_order_identifiable"])


class AlphaSelectionTest(unittest.TestCase):
    def setUp(self) -> None:
        import pandas as pd

        self.records = pd.DataFrame(
            {
                "case_id": ["a", "a", "a", "b", "b", "b"],
                "alpha": [0.5, 1.0, 2.0, 0.5, 1.0, 2.0],
                "tw_mwb_auc": [0.8, 0.8, 0.7, 0.2, 0.5, 0.4],
            }
        )

    def test_per_case_selection_uses_mwb_and_smaller_alpha_tie_break(self) -> None:
        selected = select_best_rows(self.records, (0.5, 1.0, 2.0))
        actual = dict(zip(selected["case_id"], selected["alpha"], strict=True))
        self.assertEqual(actual, {"a": 0.5, "b": 1.0})

    def test_global_selection_uses_mean_synthetic_mwb(self) -> None:
        self.assertEqual(choose_global_alpha(self.records), 1.0)


class RepositoryIntegrationTest(unittest.TestCase):
    def test_real_active_component_counts(self) -> None:
        expected = {
            "naics_sic": 410,
            "cip_soc": 46,
            "openalex": 5,
            "acs_clean": 1,
            "lodes_clean": 1,
            "gtfs_clean": 1,
            "twenty_newsgroups": 1,
        }
        for stem, count in expected.items():
            matrix, _, _ = load_matrix_csv(MATRIX_DIR / f"{stem}_original_matrix.csv")
            self.assertEqual(support_component_count(matrix), count, stem)

    def test_all_300_synthetic_active_support_graphs_are_connected(self) -> None:
        for family in FAMILIES:
            for size in SIZES:
                for seed in range(20):
                    matrix = build_case(family, size, seed)["observed"]
                    self.assertEqual(
                        support_component_count(matrix),
                        1,
                        f"{family.name}/{size.name}/seed={seed}",
                    )


if __name__ == "__main__":
    unittest.main()
