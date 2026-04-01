"""Core KS statistic and scoring functions for the CMAP algorithm.

This module implements a modified signed Kolmogorov-Smirnov statistic as described
in the original Broad Institute CMAP paper and ported from the reference R implementation
(legacy/cmap.R).  The statistic differs from scipy's ``ks_2samp`` in that it returns a
*signed* value: positive when the empirical CDF lies predominantly above the uniform
reference, negative when it lies below.

All public functions are fully vectorised over the *instance* dimension so that a
rank matrix with thousands of columns can be scored in a single NumPy call.
"""

from __future__ import annotations

import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Scalar KS statistic
# ---------------------------------------------------------------------------


def ks_statistic(v: np.ndarray, n: int) -> float:
    """Compute the signed KS statistic for a single sorted rank vector.

    This is a direct port of ``detailedMaxValue`` / ``permutedMaxValue`` /
    ``pMaxValue`` from ``legacy/cmap.R`` (all three are identical).

    The statistic measures the maximum deviation between the empirical CDF of
    the supplied rank positions and a uniform CDF over ``1..n``.

    Parameters
    ----------
    v:
        1-D array of *sorted* 1-indexed rank positions (integers).  Must be
        non-empty.
    n:
        Total number of items (genes or instances) in the rank space.

    Returns
    -------
    float
        Signed KS score.  Positive when the upward deviation ``a`` dominates
        (gene set is enriched at the top of the ranked list), negative when the
        downward deviation ``b`` dominates.

    Examples
    --------
    >>> import numpy as np
    >>> ks_statistic(np.array([1, 2, 3]), 10)
    0.6666666666666667
    """
    t = len(v)
    j = np.arange(1, t + 1, dtype=np.float64)
    v_f = v.astype(np.float64)

    a = np.max(j / t - v_f / n)
    b = np.max(v_f / n - (j - 1.0) / t)

    return float(a) if a > b else float(-b)


# ---------------------------------------------------------------------------
# Vectorised KS statistic (batch over instances)
# ---------------------------------------------------------------------------


def ks_statistic_batch(sorted_ranks: np.ndarray, n: int) -> np.ndarray:
    """Vectorised signed KS statistic over multiple instances simultaneously.

    Parameters
    ----------
    sorted_ranks:
        2-D integer array of shape ``(t, num_instances)``.  Each column
        contains the *sorted* 1-indexed rank positions for one instance.
        All columns must have the same length ``t`` (the query gene-set size).
    n:
        Total number of genes (rows) in the rank matrix.

    Returns
    -------
    np.ndarray
        1-D float array of shape ``(num_instances,)`` containing the signed KS
        score for each instance.  Returns an all-zeros array when ``t == 0``.

    Notes
    -----
    Broadcasting layout::

        j        : (t, 1)   — row indices 1..t
        ranks    : (t, K)   — sorted rank values per instance
        a        : (K,)     — max positive deviation per instance
        b        : (K,)     — max negative deviation per instance
    """
    if sorted_ranks.ndim == 1:
        sorted_ranks = sorted_ranks[:, np.newaxis]

    t, num_instances = sorted_ranks.shape

    if t == 0:
        return np.zeros(num_instances, dtype=np.float64)

    j = np.arange(1, t + 1, dtype=np.float64)[:, np.newaxis]  # (t, 1)
    ranks_f = sorted_ranks.astype(np.float64)                   # (t, K)

    a = np.max(j / t - ranks_f / n, axis=0)          # (K,)
    b = np.max(ranks_f / n - (j - 1.0) / t, axis=0)  # (K,)

    # Signed selection: use a when a > b, else -b
    return np.where(a > b, a, -b)


# ---------------------------------------------------------------------------
# Combined CMAP score
# ---------------------------------------------------------------------------


def compute_scores_batch(k_up: np.ndarray, k_down: np.ndarray) -> np.ndarray:
    """Compute the combined CMAP score from up- and down-regulation KS statistics.

    R logic (``cmap.detResults``)::

        if (k_up >= 0 && k_down >= 0 || k_up < 0 && k_down < 0):
            s_i <- 0
        else:
            s_i <- k_up - k_down

    When both KS scores have the same sign the gene set does not show a
    consistent connectivity pattern, so the combined score is set to zero.

    Parameters
    ----------
    k_up:
        1-D array of KS scores for the up-regulated gene set.
    k_down:
        1-D array of KS scores for the down-regulated gene set.

    Returns
    -------
    np.ndarray
        1-D array of combined connectivity scores, same length as inputs.
    """
    same_sign = (k_up >= 0) == (k_down >= 0)
    return np.where(same_sign, 0.0, k_up - k_down)


# ---------------------------------------------------------------------------
# Score normalisation
# ---------------------------------------------------------------------------


def normalize_scores(scores: np.ndarray) -> np.ndarray:
    """Normalise combined CMAP scores to the interval ``[-1, 1]``.

    R logic::

        p <- max(scores)
        q <- min(scores)
        for each s:
            if s >= 0: s <- s / p
            else:      s <- -(s / q)

    Edge cases handled:
    - All scores are zero → return all zeros.
    - No positive scores → positive normalisation is skipped (zeros stay zero).
    - No negative scores → negative normalisation is skipped.

    Parameters
    ----------
    scores:
        1-D array of raw combined scores from :func:`compute_scores_batch`.

    Returns
    -------
    np.ndarray
        1-D float array with values in ``[-1, 1]``.
    """
    scores = scores.astype(np.float64, copy=True)
    p = float(np.max(scores))
    q = float(np.min(scores))

    if p == 0.0 and q == 0.0:
        return scores

    result = np.empty_like(scores)

    pos_mask = scores >= 0.0
    neg_mask = ~pos_mask

    if p != 0.0:
        result[pos_mask] = scores[pos_mask] / p
    else:
        result[pos_mask] = 0.0

    if q != 0.0:
        result[neg_mask] = -(scores[neg_mask] / q)
    else:
        result[neg_mask] = 0.0

    return result


# ---------------------------------------------------------------------------
# Permutation null distribution
# ---------------------------------------------------------------------------


def ks_permutation_batch(
    n_total: int,
    sample_size: int,
    n_permutations: int = 100_000,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Generate a null distribution of absolute KS scores via random permutations.

    This is the vectorised equivalent of the ``pValue`` function in
    ``legacy/cmap.R``, which draws 100 000 random samples of ``sample_size``
    items from ``1..n_total``, sorts each sample, computes the KS statistic and
    returns the absolute values.

    Memory guard: when the full permutation matrix would exceed ~500 MB the
    work is split into chunks that each stay below the limit.

    Parameters
    ----------
    n_total:
        Total number of items in the rank space (``instances`` in R).
    sample_size:
        Number of items to draw per permutation (``lens`` in R).
    n_permutations:
        Number of random permutations.  Defaults to 100 000 (R default).
    rng:
        Optional seeded :class:`numpy.random.Generator`.  When ``None`` a
        default generator is created via :func:`numpy.random.default_rng`.

    Returns
    -------
    np.ndarray
        1-D float array of shape ``(n_permutations,)`` containing absolute KS
        scores.

    Notes
    -----
    Each column of the permutation matrix is a sorted random sample drawn
    *without* replacement (matching ``sample.int`` in R which defaults to
    ``replace=FALSE``).
    """
    if rng is None:
        rng = np.random.default_rng()

    # Bytes per permutation column: sample_size * 8 (float64 after sort)
    bytes_per_col = sample_size * 8
    # Target: stay under 500 MB per chunk
    chunk_size = max(1, 500_000_000 // bytes_per_col)

    results: list[np.ndarray] = []
    remaining = n_permutations

    while remaining > 0:
        batch = min(chunk_size, remaining)

        # Draw `batch` columns of `sample_size` random integers from 1..n_total
        # without replacement.  shape: (sample_size, batch)
        samples = np.stack(
            [rng.choice(n_total, size=sample_size, replace=False) + 1 for _ in range(batch)],
            axis=1,
        ).astype(np.int32)

        # Sort each column in-place
        samples.sort(axis=0)

        ks_vals = ks_statistic_batch(samples, n_total)  # (batch,)
        results.append(np.abs(ks_vals))
        remaining -= batch

    return np.concatenate(results)


# ---------------------------------------------------------------------------
# Non-null percentage
# ---------------------------------------------------------------------------


def non_null_pct(scores: np.ndarray) -> float:
    """Compute the non-null percentage of a score vector.

    R logic (``nonNull`` in ``legacy/cmap.R``)::

        lessNull    <- count(scores < 0)
        greaterNull <- count(scores > 0)
        if both 0: return 0
        return max(lessNull, greaterNull) * 100 / length(scores)

    Parameters
    ----------
    scores:
        1-D array of CMAP scores (typically the ``Score`` column of detailed
        results for a single drug name).

    Returns
    -------
    float
        Percentage of scores that are non-zero, weighted toward whichever tail
        (positive or negative) is larger.  Returns ``0.0`` when all scores are
        exactly zero.
    """
    less_null = int(np.sum(scores < 0))
    greater_null = int(np.sum(scores > 0))

    if less_null == 0 and greater_null == 0:
        return 0.0

    return max(less_null, greater_null) * 100.0 / len(scores)
