"""Permuted results computation for the CMAP algorithm.

This module implements the equivalent of ``cmap.permResults`` from
``legacy/cmap.R``.  It computes enrichment, p-value, specificity and
non-null percentage for each unique drug name (``"by_name"``) or each
unique drug-name + cell-line combination (``"by_cell"``).

The R code groups instances by name, finds their row positions within the
sorted detailed results, and computes a KS statistic on those positions.
This Python port preserves the exact same mathematical logic while
replacing the R ``for`` loop with partially vectorised group operations.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from pycmap._types import CmapData, DetailedResults, PermutedResults
from pycmap.core import ks_statistic
from pycmap.core import non_null_pct as _non_null_pct

# ---------------------------------------------------------------------------
# Per-group helper functions
# ---------------------------------------------------------------------------


def compute_p_value(
    enrichment: float,
    n: int,
    mean_score: float,
    non_null: float,
    p_value_dist: dict[str, np.ndarray],
) -> float:
    """Compute the empirical p-value from a pre-computed null distribution.

    Returns the sentinel value ``100.0`` when the instance group is too
    small, has a zero mean, or has a non-null percentage ≤ 50 — exactly
    matching the R implementation's behaviour.

    Parameters
    ----------
    enrichment:
        Observed enrichment score for this group.
    n:
        Number of instances in the group (frequency).
    mean_score:
        Mean of the detailed scores for this group's instances.
    non_null:
        Non-null percentage for this group.
    p_value_dist:
        Pre-computed null distributions keyed by frequency string (e.g.
        ``"2"``, ``"15"``).  Each value is a 1-D array of absolute KS
        scores from random permutations.

    Returns
    -------
    float
        Empirical p-value, or ``100.0`` as a sentinel.
    """
    if n == 1 or mean_score == 0.0 or non_null <= 50.0:
        return 100.0

    key = str(n)
    if key not in p_value_dist:
        return 100.0

    null_dist = p_value_dist[key]
    return float(np.sum(null_dist >= abs(enrichment)) / len(null_dist))


def compute_specificity(enrichment: float, msig_scores: np.ndarray) -> float:
    """Compute specificity as a percentile rank within same-direction MSigDB scores.

    The specificity measures how specific the enrichment is to *this* drug
    compared to all MSigDB gene-set queries.  Lower values indicate higher
    specificity (more extreme than other pathways).

    Parameters
    ----------
    enrichment:
        Observed enrichment score for this group.
    msig_scores:
        1-D array of enrichment scores for this instance across all MSigDB
        gene sets (one row of the pre-computed specificity matrix).

    Returns
    -------
    float
        Specificity in ``[0, 1]``.

    Notes
    -----
    R logic from ``legacy/cmap.R`` lines 264-296::

        if enrichment == 0: specificity = 0
        elif enrichment < 0:
            negatives = msig_scores[msig_scores < 0]
            specificity = count(enrichment >= negatives) / len(negatives)
        elif enrichment > 0:
            positives = msig_scores[msig_scores > 0]
            specificity = 1 - count(enrichment >= positives) / len(positives)
    """
    if enrichment == 0.0:
        return 0.0

    if enrichment < 0.0:
        neg_scores = msig_scores[msig_scores < 0]
        if len(neg_scores) == 0:
            return 0.0
        return float(np.sum(enrichment >= neg_scores) / len(neg_scores))

    # enrichment > 0
    pos_scores = msig_scores[msig_scores > 0]
    if len(pos_scores) == 0:
        return 0.0
    return float(1.0 - np.sum(enrichment >= pos_scores) / len(pos_scores))


# ---------------------------------------------------------------------------
# Main permuted results computation
# ---------------------------------------------------------------------------


def _add_cell_line_labels(
    detailed: DetailedResults,
    data: CmapData,
) -> tuple[DetailedResults, list[str]]:
    """Re-label detailed results with 'name - cellline' and re-sort by instance ID.

    This matches the R ``addCellLine`` function which:
    1. Sorts by instance ID
    2. Assigns new row names as 'instance_name - cell_line'
    3. Re-sorts by score descending, then ks_up descending
    """
    # Build a lookup from instance_id -> cell_line
    id_to_cell: dict[int, str] = {}
    id_to_name: dict[int, str] = {}
    for i, inst_id in enumerate(data.instance_ids):
        id_to_cell[int(inst_id)] = data.cell_lines[i]
        id_to_name[int(inst_id)] = data.instance_names[i]

    # Create "name - cellline" labels for each row of detailed results
    new_names = [
        f"{id_to_name.get(int(iid), '?')} - {id_to_cell.get(int(iid), '?')}"
        for iid in detailed.instance_ids
    ]

    # Re-sort by score desc, ks_up desc (same as original)
    order = np.lexsort((-detailed.ks_up, -detailed.scores))

    relabeled = DetailedResults(
        instance_ids=detailed.instance_ids[order],
        scores=detailed.scores[order],
        ks_up=detailed.ks_up[order],
        ks_down=detailed.ks_down[order],
        instance_names=[new_names[i] for i in order],
    )

    return relabeled, data.unique_name_cell


def compute_permuted_results(
    detailed: DetailedResults,
    data: CmapData,
    mode: Literal["by_name", "by_cell"] = "by_name",
    p_value_dist: dict[str, np.ndarray] | None = None,
    specificity_matrix: np.ndarray | None = None,
) -> PermutedResults:
    """Compute permuted CMAP results aggregated by drug name or name+cell line.

    This is the equivalent of ``cmap.permResults`` from ``legacy/cmap.R``.
    For each unique group (drug name or drug-cell combination):

    1. Find all row positions matching that name in the sorted detailed results.
    2. Compute **enrichment** = KS statistic on those row positions.
    3. Compute **mean** of detailed scores for matching instances.
    4. Compute **non-null %** = percentage of the larger sign-group.
    5. Compute **p-value** from pre-computed null distribution.
    6. Compute **specificity** from MSigDB cross-enrichment matrix.

    Parameters
    ----------
    detailed:
        Output from :func:`~pycmap.detailed.compute_detailed_results`.
    data:
        Loaded CMAP data.
    mode:
        ``"by_name"`` groups by drug name; ``"by_cell"`` groups by
        drug name + cell line.
    p_value_dist:
        Pre-computed null distributions keyed by frequency string.
        If ``None``, all p-values will be set to the sentinel ``100.0``.
    specificity_matrix:
        Pre-computed specificity matrix of shape
        ``(n_unique_groups, n_msigdb_sets)``.  Row order must match
        ``data.unique_names`` (for ``by_name``) or
        ``data.unique_name_cell`` (for ``by_cell``).
        If ``None``, all specificities will be ``0.0``.

    Returns
    -------
    PermutedResults
        Sorted by p-value ascending, then enrichment descending (matching R).
    """
    if p_value_dist is None:
        p_value_dist = {}

    # --- determine grouping ---------------------------------------------------
    if mode == "by_cell":
        det, unique_names = _add_cell_line_labels(detailed, data)
    else:
        det = detailed
        unique_names = data.unique_names

    n_instances = len(det.scores)

    # Build a name -> list of row indices mapping
    name_to_rows: dict[str, list[int]] = {}
    for i, name in enumerate(det.instance_names):
        name_to_rows.setdefault(name, []).append(i)

    # --- compute per-group metrics -------------------------------------------
    n_groups = len(unique_names)
    means = np.empty(n_groups, dtype=np.float64)
    ns = np.empty(n_groups, dtype=np.int64)
    enrichments = np.empty(n_groups, dtype=np.float64)
    p_values = np.empty(n_groups, dtype=np.float64)
    specificities = np.empty(n_groups, dtype=np.float64)
    non_null_pcts = np.empty(n_groups, dtype=np.float64)

    for i, name in enumerate(unique_names):
        rows = name_to_rows.get(name, [])
        n = len(rows)
        ns[i] = n

        if n == 0:
            means[i] = 0.0
            enrichments[i] = 0.0
            p_values[i] = 100.0
            specificities[i] = 0.0
            non_null_pcts[i] = 0.0
            continue

        # Row positions are 1-indexed for KS (matching R's which() returning 1-indexed)
        row_positions = np.array(sorted(rows), dtype=np.float64) + 1.0
        instance_scores = det.scores[rows]

        enrichments[i] = ks_statistic(row_positions, n_instances)
        means[i] = float(np.mean(instance_scores))
        non_null_pcts[i] = _non_null_pct(instance_scores)

        # P-value
        p_values[i] = compute_p_value(
            enrichments[i], n, means[i], non_null_pcts[i], p_value_dist
        )

        # Specificity
        if specificity_matrix is not None and i < specificity_matrix.shape[0]:
            specificities[i] = compute_specificity(
                enrichments[i], specificity_matrix[i, :]
            )
        else:
            specificities[i] = 0.0

    # --- sort: p-value ascending, then enrichment descending -----------------
    # np.lexsort sorts by last key first
    order = np.lexsort((-enrichments, p_values))

    return PermutedResults(
        mean=means[order],
        n=ns[order],
        enrichment=enrichments[order],
        p_value=p_values[order],
        specificity=specificities[order],
        non_null_pct=non_null_pcts[order],
        names=[unique_names[i] for i in order],
    )
