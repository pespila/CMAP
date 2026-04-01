"""Detailed results computation for the CMAP algorithm.

This module implements the equivalent of ``cmap.detResults`` from
``legacy/cmap.R``.  For each drug instance (column in the rank matrix) it
computes a signed KS statistic for the up- and down-regulated gene sets,
combines them into a single connectivity score, and normalises the scores
to the interval ``[-1, 1]``.

The entire computation is **fully vectorised** over instances — there is no
Python-level loop over the ~6 100 columns of the rank matrix.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

from pycmap._types import CmapData, DetailedResults
from pycmap.core import compute_scores_batch, ks_statistic_batch, normalize_scores
from pycmap.data import build_gene_index, read_grp


def _resolve_genes(
    genes: list[str] | str | Path,
) -> list[str]:
    """Accept a gene list, a .grp file path, or a Path and return a list of gene IDs."""
    if isinstance(genes, (str, Path)):
        return read_grp(genes)
    return list(genes)


def _lookup_gene_indices(
    genes: list[str],
    gene_index: dict[str, int],
    label: str,
) -> list[int]:
    """Map gene names to row indices, warning about any missing genes."""
    found: list[int] = []
    missing: list[str] = []
    for g in genes:
        idx = gene_index.get(g)
        if idx is not None:
            found.append(idx)
        else:
            missing.append(g)

    if missing:
        warnings.warn(
            f"{len(missing)} {label} gene(s) not found in rank matrix and will be "
            f"ignored: {missing[:5]}{'...' if len(missing) > 5 else ''}",
            stacklevel=3,
        )

    if not found:
        raise ValueError(
            f"None of the {label} genes were found in the rank matrix. "
            "Check that gene identifiers match the rank matrix row names."
        )

    return found


def compute_detailed_results(
    up_genes: list[str] | str | Path,
    down_genes: list[str] | str | Path,
    data: CmapData,
) -> DetailedResults:
    """Compute detailed CMAP results for all instances.

    This is the vectorised equivalent of ``cmap.detResults`` in
    ``legacy/cmap.R``.  The entire rank matrix is scored in two batch KS
    calls (one for up-regulated genes, one for down-regulated genes)
    followed by vectorised score combination and normalisation.

    Parameters
    ----------
    up_genes:
        Up-regulated gene identifiers — a ``list[str]``, a path to a
        ``.grp`` file, or a :class:`~pathlib.Path`.
    down_genes:
        Down-regulated gene identifiers (same input types as *up_genes*).
    data:
        Loaded CMAP data (see :func:`~pycmap.data.load_cmap_data`).

    Returns
    -------
    DetailedResults
        Sorted by ``scores`` descending, then ``ks_up`` descending for ties
        (matching the R implementation's ``order`` call).

    Raises
    ------
    ValueError
        When none of the supplied genes are found in the rank matrix.

    Warnings
    --------
    UserWarning
        When some (but not all) genes are missing from the rank matrix.

    Examples
    --------
    >>> from pycmap.data import load_cmap_data
    >>> data = load_cmap_data("/path/to/cmap/data")
    >>> result = compute_detailed_results("readUp.grp", "readDown.grp", data)
    >>> df = result.to_dataframe()
    >>> df.head()
    """
    # --- resolve gene lists --------------------------------------------------
    up_list = _resolve_genes(up_genes)
    down_list = _resolve_genes(down_genes)

    gene_index = build_gene_index(data.gene_names)

    up_indices = _lookup_gene_indices(up_list, gene_index, "up-regulated")
    down_indices = _lookup_gene_indices(down_list, gene_index, "down-regulated")

    n_genes = data.rank_matrix.shape[0]

    # --- extract and sort sub-matrices ---------------------------------------
    # Shape: (t_up, n_instances) and (t_down, n_instances)
    up_ranks = np.sort(data.rank_matrix[up_indices, :], axis=0)
    down_ranks = np.sort(data.rank_matrix[down_indices, :], axis=0)

    # --- vectorised KS -------------------------------------------------------
    k_up = ks_statistic_batch(up_ranks, n_genes)
    k_down = ks_statistic_batch(down_ranks, n_genes)

    # --- combined score and normalisation ------------------------------------
    raw_scores = compute_scores_batch(k_up, k_down)
    scores = normalize_scores(raw_scores)

    # --- sort: score descending, then k_up descending for ties ---------------
    # np.lexsort sorts by last key first, so order is (-k_up, -scores)
    order = np.lexsort((-k_up, -scores))

    return DetailedResults(
        instance_ids=data.instance_ids[order],
        scores=scores[order],
        ks_up=k_up[order],
        ks_down=k_down[order],
        instance_names=[data.instance_names[i] for i in order],
    )
