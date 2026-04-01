"""Tests for pycmap.permuted — compute_permuted_results, compute_p_value,
compute_specificity, and non_null_pct integration.
"""

from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose

from pycmap._types import CmapData, PermutedResults
from pycmap.core import non_null_pct
from pycmap.detailed import compute_detailed_results
from pycmap.permuted import compute_p_value, compute_permuted_results, compute_specificity

# ---------------------------------------------------------------------------
# non_null_pct (via pycmap.core) — additional edge cases
# ---------------------------------------------------------------------------


def test_non_null_pct_single_positive() -> None:
    """A length-1 array with one positive value → 100 %."""
    assert_allclose(non_null_pct(np.array([0.5])), 100.0)


def test_non_null_pct_single_negative() -> None:
    """A length-1 array with one negative value → 100 %."""
    assert_allclose(non_null_pct(np.array([-0.5])), 100.0)


def test_non_null_pct_single_zero() -> None:
    """A length-1 array with zero → 0 %."""
    assert non_null_pct(np.array([0.0])) == 0.0


def test_non_null_pct_equal_split() -> None:
    """Equal number of positive and negative (and some zeros) → max / total."""
    scores = np.array([0.1, -0.1, 0.2, -0.2, 0.0])
    # 2 positive, 2 negative → max=2 → 2/5*100=40.0
    assert_allclose(non_null_pct(scores), 40.0)


# ---------------------------------------------------------------------------
# compute_p_value — sentinel conditions
# ---------------------------------------------------------------------------


def test_compute_p_value_sentinel_n1() -> None:
    """n=1 must return the sentinel 100.0 regardless of other args."""
    dist = {"1": np.array([0.1, 0.2, 0.3])}
    result = compute_p_value(
        enrichment=0.5, n=1, mean_score=0.4, non_null=80.0, p_value_dist=dist
    )
    assert result == 100.0


def test_compute_p_value_sentinel_mean_zero() -> None:
    """mean_score=0 must return the sentinel 100.0."""
    dist = {"3": np.array([0.1, 0.2, 0.3])}
    result = compute_p_value(
        enrichment=0.5, n=3, mean_score=0.0, non_null=80.0, p_value_dist=dist
    )
    assert result == 100.0


def test_compute_p_value_sentinel_non_null_low() -> None:
    """non_null <= 50 must return the sentinel 100.0."""
    dist = {"4": np.array([0.1, 0.2, 0.3, 0.4])}
    result = compute_p_value(
        enrichment=0.5, n=4, mean_score=0.3, non_null=50.0, p_value_dist=dist
    )
    assert result == 100.0


def test_compute_p_value_sentinel_non_null_exactly_50() -> None:
    """non_null exactly 50 still triggers the sentinel (condition is <= 50)."""
    dist = {"5": np.array([0.1, 0.2, 0.3])}
    result = compute_p_value(
        enrichment=0.3, n=5, mean_score=0.2, non_null=50.0, p_value_dist=dist
    )
    assert result == 100.0


def test_compute_p_value_sentinel_missing_key() -> None:
    """Missing key in p_value_dist must return the sentinel 100.0."""
    dist: dict[str, np.ndarray] = {}
    result = compute_p_value(
        enrichment=0.5, n=4, mean_score=0.3, non_null=80.0, p_value_dist=dist
    )
    assert result == 100.0


def test_compute_p_value_normal() -> None:
    """With a concrete null distribution the p-value must equal the empirical tail proportion."""
    # null distribution: 10 values uniformly in [0, 1]
    null = np.linspace(0.0, 1.0, 11)  # [0.0, 0.1, ..., 1.0]
    dist = {"4": null}

    enrichment = 0.55  # abs(0.55) = 0.55
    # Proportion of null >= 0.55: values are 0.6, 0.7, 0.8, 0.9, 1.0 → 5 out of 11
    expected_p = float(np.sum(null >= abs(enrichment)) / len(null))

    result = compute_p_value(
        enrichment=enrichment, n=4, mean_score=0.4, non_null=80.0, p_value_dist=dist
    )
    assert_allclose(result, expected_p, rtol=1e-12)


def test_compute_p_value_zero_enrichment_included() -> None:
    """Enrichment=0 with positive mean and non_null > 50 returns an empirical p-value."""
    null = np.array([0.0, 0.1, 0.2, 0.3])
    dist = {"3": null}
    # abs(0.0) = 0 → all 4 values in null >= 0 → p = 4/4 = 1.0
    result = compute_p_value(
        enrichment=0.0, n=3, mean_score=0.3, non_null=80.0, p_value_dist=dist
    )
    assert_allclose(result, 1.0)


# ---------------------------------------------------------------------------
# compute_specificity
# ---------------------------------------------------------------------------


def test_compute_specificity_zero_enrichment() -> None:
    """enrichment=0 must always return 0.0."""
    msig = np.array([-0.5, -0.2, 0.1, 0.4, 0.8])
    assert compute_specificity(0.0, msig) == 0.0


def test_compute_specificity_negative() -> None:
    """Verify computation for a negative enrichment value."""
    msig = np.array([-0.5, -0.3, -0.1, 0.2, 0.4, 0.6])
    # neg_scores = [-0.5, -0.3, -0.1]
    # enrichment=-0.2 >= [-0.5,-0.3,-0.1]: [T, T, F] → 2/3
    result = compute_specificity(-0.2, msig)
    assert_allclose(result, 2.0 / 3.0, rtol=1e-12)


def test_compute_specificity_positive() -> None:
    """Verify computation for a positive enrichment value."""
    msig = np.array([-0.5, -0.3, -0.1, 0.2, 0.4, 0.6])
    # pos_scores = [0.2, 0.4, 0.6]
    # enrichment=0.5 >= [0.2, 0.4, 0.6]: [T, T, F] → 2/3 → specificity = 1 - 2/3 = 1/3
    result = compute_specificity(0.5, msig)
    assert_allclose(result, 1.0 / 3.0, rtol=1e-12)


def test_compute_specificity_no_same_direction_scores() -> None:
    """When there are no same-direction MSigDB scores the specificity is 0.0."""
    msig_all_positive = np.array([0.1, 0.3, 0.5])
    # enrichment < 0 but no negatives in msig
    assert compute_specificity(-0.4, msig_all_positive) == 0.0

    msig_all_negative = np.array([-0.1, -0.3, -0.5])
    # enrichment > 0 but no positives in msig
    assert compute_specificity(0.4, msig_all_negative) == 0.0


def test_compute_specificity_highest_enrichment_is_zero() -> None:
    """Enrichment that beats all same-direction scores → specificity == 0.0 (positive)."""
    msig = np.array([-0.5, 0.1, 0.2, 0.3])
    # enrichment=0.9 > all pos_scores → np.sum(0.9 >= [0.1,0.2,0.3])=3 → 1-3/3=0.0
    result = compute_specificity(0.9, msig)
    assert_allclose(result, 0.0, atol=1e-12)


def test_compute_specificity_lowest_negative_enrichment() -> None:
    """Enrichment that beats none of the neg scores → specificity == 0.0 (negative)."""
    msig = np.array([-0.5, -0.3, -0.1, 0.2])
    # enrichment=-0.6 < all neg_scores → np.sum(-0.6 >= [-0.5,-0.3,-0.1])=0 → 0/3=0.0
    result = compute_specificity(-0.6, msig)
    assert_allclose(result, 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# compute_permuted_results — shape and ordering
# ---------------------------------------------------------------------------


def test_permuted_results_shape(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Output PermutedResults arrays must have length == number of unique names."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data, mode="by_name")

    n_unique = len(cmap_data.unique_names)
    assert isinstance(permuted, PermutedResults)
    assert len(permuted.mean) == n_unique
    assert len(permuted.n) == n_unique
    assert len(permuted.enrichment) == n_unique
    assert len(permuted.p_value) == n_unique
    assert len(permuted.specificity) == n_unique
    assert len(permuted.non_null_pct) == n_unique
    assert len(permuted.names) == n_unique


def test_permuted_results_by_cell_shape(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """by_cell mode must produce length == number of unique name-cell combinations."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data, mode="by_cell")

    n_unique_cell = len(cmap_data.unique_name_cell)
    assert len(permuted.names) == n_unique_cell


def test_permuted_results_sorted_by_p(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """p_value array must be in non-decreasing (ascending) order."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data, mode="by_name")
    assert np.all(np.diff(permuted.p_value) >= 0), (
        "p_values are not in ascending order"
    )


def test_permuted_results_n_counts_correct(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """n counts must sum to the total number of instances."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data, mode="by_name")
    assert int(np.sum(permuted.n)) == cmap_data.rank_matrix.shape[1]


def test_permuted_results_sentinel_p_values_no_dist(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """With no p_value_dist supplied, all p-values must be the sentinel 100.0."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data, mode="by_name", p_value_dist=None)
    assert_allclose(permuted.p_value, 100.0)


def test_permuted_results_specificity_zero_no_matrix(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """With no specificity_matrix supplied, all specificities must be 0.0."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(
        detailed, cmap_data, mode="by_name", specificity_matrix=None
    )
    assert_allclose(permuted.specificity, 0.0)


def test_permuted_results_to_dataframe_columns(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """to_dataframe() must produce exactly the expected column names."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    permuted = compute_permuted_results(detailed, cmap_data)
    df = permuted.to_dataframe()
    assert list(df.columns) == ["mean", "n", "enrichment", "p", "specificity", "% non-null"]


def test_permuted_results_deterministic(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Two consecutive calls with the same inputs must produce identical results."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    p1 = compute_permuted_results(detailed, cmap_data, mode="by_name")
    p2 = compute_permuted_results(detailed, cmap_data, mode="by_name")
    assert_allclose(p1.mean, p2.mean)
    assert_allclose(p1.enrichment, p2.enrichment)
    assert list(p1.names) == list(p2.names)


def test_permuted_results_with_p_value_dist(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """With a non-empty p_value_dist, groups with n>1 that qualify may get real p-values."""
    detailed = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)

    # Build a null distribution for groups of size 4 (each drug appears 4 times)
    null = np.abs(
        np.array([
            float(np.random.default_rng(i).uniform(-1, 1)) for i in range(1000)
        ])
    )
    p_value_dist = {"4": null}

    permuted = compute_permuted_results(
        detailed, cmap_data, mode="by_name", p_value_dist=p_value_dist
    )
    # Result must have the right length — values may still be 100.0 if mean==0 or non_null<=50
    assert len(permuted.p_value) == len(cmap_data.unique_names)
    assert np.all(permuted.p_value >= 0.0)
    assert np.all(permuted.p_value <= 100.0)
