"""Tests for pycmap.core — KS statistics, score combination, normalisation,
permutation null distribution, and non-null percentage.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycmap.core import (
    compute_scores_batch,
    ks_permutation_batch,
    ks_statistic,
    ks_statistic_batch,
    non_null_pct,
    normalize_scores,
)

# ---------------------------------------------------------------------------
# ks_statistic — scalar
# ---------------------------------------------------------------------------


def test_ks_statistic_known_values() -> None:
    """v=[1,2,3], n=10 should return 0.7 (upward deviation dominates)."""
    v = np.array([1, 2, 3])
    result = ks_statistic(v, 10)
    # a = max(j/t - v/n) = max([1/3-0.1, 2/3-0.2, 1-0.3]) = max([0.233,0.467,0.7]) = 0.7
    # b = max(v/n - (j-1)/t) = max([0.1, -0.133, -0.367]) = 0.1
    # a > b → return a = 0.7
    assert_allclose(result, 0.7, rtol=1e-10)


def test_ks_statistic_negative_result() -> None:
    """v=[8,9,10], n=10 should return a negative score (genes at bottom)."""
    v = np.array([8, 9, 10])
    result = ks_statistic(v, 10)
    # a = max(j/3 - v/10) = max([-0.467,-0.233,0.0]) = 0.0
    # b = max(v/10 - (j-1)/3) = max([0.8, 0.567, 0.333]) = 0.8
    # a <= b → return -b = -0.8
    assert_allclose(result, -0.8, rtol=1e-10)


def test_ks_statistic_single_element_top() -> None:
    """Single element at position 1 in n=10 — strong upward enrichment."""
    v = np.array([1])
    result = ks_statistic(v, 10)
    # a = max(1/1 - 1/10) = 0.9
    # b = max(1/10 - 0/1) = 0.1
    # a > b → 0.9
    assert_allclose(result, 0.9, rtol=1e-10)


def test_ks_statistic_single_element_bottom() -> None:
    """Single element at last position n in n — strong downward enrichment."""
    n = 10
    v = np.array([n])
    result = ks_statistic(v, n)
    # a = max(1/1 - 10/10) = 0.0
    # b = max(10/10 - 0/1) = 1.0
    # a <= b → -1.0
    assert_allclose(result, -1.0, rtol=1e-10)


# ---------------------------------------------------------------------------
# ks_statistic_batch — vectorised
# ---------------------------------------------------------------------------


def test_ks_statistic_batch_matches_scalar() -> None:
    """Batch result must equal a loop of scalar calls for every column."""
    rng = np.random.default_rng(7)
    n = 50
    t = 8
    k = 15
    # Build a (t, K) matrix of sorted rank draws
    cols = [np.sort(rng.choice(n, size=t, replace=False) + 1) for _ in range(k)]
    sorted_ranks = np.stack(cols, axis=1).astype(np.int32)

    batch_result = ks_statistic_batch(sorted_ranks, n)
    scalar_results = np.array([ks_statistic(sorted_ranks[:, j], n) for j in range(k)])

    assert_allclose(batch_result, scalar_results, rtol=1e-12)


def test_ks_statistic_batch_1d_input() -> None:
    """A 1-D input should be treated as a single-column matrix."""
    v = np.array([1, 2, 3])
    scalar = ks_statistic(v, 10)
    batch = ks_statistic_batch(v, 10)
    assert batch.shape == (1,)
    assert_allclose(batch[0], scalar, rtol=1e-12)


def test_ks_statistic_batch_empty_returns_zeros() -> None:
    """Empty (t=0) sorted_ranks should return an all-zeros array."""
    empty = np.empty((0, 5), dtype=np.int32)
    result = ks_statistic_batch(empty, 100)
    assert result.shape == (5,)
    assert_allclose(result, np.zeros(5))


def test_ks_statistic_batch_shape() -> None:
    """Output shape must equal (num_instances,)."""
    rng = np.random.default_rng(99)
    n, t, k = 100, 10, 20
    cols = [np.sort(rng.choice(n, size=t, replace=False) + 1) for _ in range(k)]
    sorted_ranks = np.stack(cols, axis=1).astype(np.int32)
    result = ks_statistic_batch(sorted_ranks, n)
    assert result.shape == (k,)


# ---------------------------------------------------------------------------
# compute_scores_batch
# ---------------------------------------------------------------------------


def test_compute_scores_same_sign_zero_positive() -> None:
    """k_up=0.5, k_down=0.3 — both positive → combined score must be 0."""
    k_up = np.array([0.5])
    k_down = np.array([0.3])
    result = compute_scores_batch(k_up, k_down)
    assert_allclose(result, np.array([0.0]))


def test_compute_scores_same_sign_zero_negative() -> None:
    """k_up=-0.4, k_down=-0.2 — both negative → combined score must be 0."""
    k_up = np.array([-0.4])
    k_down = np.array([-0.2])
    result = compute_scores_batch(k_up, k_down)
    assert_allclose(result, np.array([0.0]))


def test_compute_scores_opposite_sign() -> None:
    """k_up=0.5, k_down=-0.3 → score = k_up - k_down = 0.5 - (-0.3) = 0.8."""
    k_up = np.array([0.5])
    k_down = np.array([-0.3])
    result = compute_scores_batch(k_up, k_down)
    assert_allclose(result, np.array([0.8]), rtol=1e-12)


def test_compute_scores_opposite_sign_negative_up() -> None:
    """k_up=-0.5, k_down=0.3 → score = -0.5 - 0.3 = -0.8."""
    k_up = np.array([-0.5])
    k_down = np.array([0.3])
    result = compute_scores_batch(k_up, k_down)
    assert_allclose(result, np.array([-0.8]), rtol=1e-12)


def test_compute_scores_batch_mixed_array() -> None:
    """Batch with a mix of same-sign and opposite-sign pairs."""
    k_up = np.array([0.5, -0.4, 0.3, -0.2])
    k_down = np.array([0.3, -0.1, -0.2, 0.6])
    result = compute_scores_batch(k_up, k_down)
    expected = np.array([0.0, 0.0, 0.5, -0.8])
    assert_allclose(result, expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# normalize_scores
# ---------------------------------------------------------------------------


def test_normalize_scores_known_values() -> None:
    """[0.8, -0.4, 0.4, -0.2] → [1.0, -1.0, 0.5, -0.5]."""
    scores = np.array([0.8, -0.4, 0.4, -0.2])
    result = normalize_scores(scores)
    assert_allclose(result, np.array([1.0, -1.0, 0.5, -0.5]), rtol=1e-12)


def test_normalize_scores_all_zeros() -> None:
    """All-zero input must return all zeros."""
    scores = np.zeros(5)
    result = normalize_scores(scores)
    assert_allclose(result, np.zeros(5))


def test_normalize_scores_only_positive() -> None:
    """Only positive scores — max maps to 1.0."""
    scores = np.array([0.2, 0.4, 0.6, 1.0])
    result = normalize_scores(scores)
    assert_allclose(result, np.array([0.2, 0.4, 0.6, 1.0]), rtol=1e-12)


def test_normalize_scores_only_negative() -> None:
    """Only negative scores — min maps to -1.0."""
    scores = np.array([-0.5, -1.0, -0.25])
    result = normalize_scores(scores)
    # q = -1.0; result = -(score / -1.0) for negatives
    assert_allclose(result, np.array([-0.5, -1.0, -0.25]), rtol=1e-12)


def test_normalize_scores_output_range() -> None:
    """All normalised scores must lie within [-1, 1]."""
    rng = np.random.default_rng(11)
    raw = rng.uniform(-5.0, 5.0, size=200)
    result = normalize_scores(raw)
    assert np.all(result >= -1.0)
    assert np.all(result <= 1.0)


# ---------------------------------------------------------------------------
# non_null_pct
# ---------------------------------------------------------------------------


def test_non_null_pct_mixed() -> None:
    """[0.3, -0.2, 0.1, 0.0, -0.5] — 2 positive, 2 negative → max=2 → 2/5*100=40.0."""
    scores = np.array([0.3, -0.2, 0.1, 0.0, -0.5])
    result = non_null_pct(scores)
    assert_allclose(result, 40.0)


def test_non_null_pct_all_zero() -> None:
    """All-zero scores must return 0.0."""
    scores = np.zeros(10)
    result = non_null_pct(scores)
    assert result == 0.0


def test_non_null_pct_all_positive() -> None:
    """All strictly positive scores → 100 %."""
    scores = np.array([0.1, 0.2, 0.5, 1.0])
    result = non_null_pct(scores)
    assert_allclose(result, 100.0)


def test_non_null_pct_all_negative() -> None:
    """All strictly negative scores → 100 %."""
    scores = np.array([-0.1, -0.5, -1.0])
    result = non_null_pct(scores)
    assert_allclose(result, 100.0)


def test_non_null_pct_more_positive() -> None:
    """More positive than negative — max is the positive count."""
    scores = np.array([0.1, 0.2, 0.3, -0.1, 0.0])
    # 3 positive, 1 negative → max=3 → 3/5*100=60.0
    result = non_null_pct(scores)
    assert_allclose(result, 60.0)


def test_non_null_pct_more_negative() -> None:
    """More negative than positive — max is the negative count."""
    scores = np.array([0.1, -0.2, -0.3, -0.4, 0.0])
    # 1 positive, 3 negative → max=3 → 3/5*100=60.0
    result = non_null_pct(scores)
    assert_allclose(result, 60.0)


# ---------------------------------------------------------------------------
# ks_permutation_batch
# ---------------------------------------------------------------------------


def test_ks_permutation_batch_shape() -> None:
    """Output shape must be (n_permutations,) and all values non-negative."""
    rng = np.random.default_rng(123)
    result = ks_permutation_batch(n_total=50, sample_size=5, n_permutations=200, rng=rng)
    assert result.shape == (200,)
    assert np.all(result >= 0.0)


def test_ks_permutation_batch_reproducible() -> None:
    """Same seed must produce identical results on two separate calls."""
    result_a = ks_permutation_batch(
        n_total=30, sample_size=4, n_permutations=100,
        rng=np.random.default_rng(77),
    )
    result_b = ks_permutation_batch(
        n_total=30, sample_size=4, n_permutations=100,
        rng=np.random.default_rng(77),
    )
    assert_allclose(result_a, result_b)


def test_ks_permutation_batch_different_seeds_differ() -> None:
    """Different seeds should (with overwhelming probability) produce different results."""
    result_a = ks_permutation_batch(
        n_total=50, sample_size=5, n_permutations=50,
        rng=np.random.default_rng(1),
    )
    result_b = ks_permutation_batch(
        n_total=50, sample_size=5, n_permutations=50,
        rng=np.random.default_rng(2),
    )
    assert not np.allclose(result_a, result_b)


@pytest.mark.parametrize("n_perms", [1, 10, 100])
def test_ks_permutation_batch_various_sizes(n_perms: int) -> None:
    """Output length equals the requested number of permutations."""
    rng = np.random.default_rng(0)
    result = ks_permutation_batch(n_total=20, sample_size=3, n_permutations=n_perms, rng=rng)
    assert len(result) == n_perms
