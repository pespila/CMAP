"""pycmap-score: Connectivity Map (CMAP) algorithm for drug-gene connectivity scoring.

This package implements the Connectivity Map algorithm described in:

    Lamb J, Crawford ED, Peck D, et al. The Connectivity Map: Using
    Gene-Expression Signatures to Connect Small Molecules, Genes, and
    Disease. *Science* 313(5795):1929-1935, 2006.
    DOI: `10.1126/science.1132939 <https://doi.org/10.1126/science.1132939>`_

Quick start
-----------
>>> from pycmap import cmap
>>> detailed, by_name, by_cell = cmap("readUp.grp", "readDown.grp", data_dir="/path/to/data")
>>> by_name.to_dataframe().head()

See :func:`cmap` for the full API.
"""

from __future__ import annotations

from pathlib import Path

from pycmap._types import CmapData, DetailedResults, PermutedResults
from pycmap.core import ks_statistic, ks_statistic_batch, non_null_pct
from pycmap.data import load_cmap_data, load_rank_matrix, read_gmt, read_grp
from pycmap.detailed import compute_detailed_results
from pycmap.install import load_p_values, load_specificity
from pycmap.permuted import compute_permuted_results

__version__ = "0.1.0"

__all__ = [
    # High-level
    "cmap",
    # Core
    "ks_statistic",
    "ks_statistic_batch",
    "non_null_pct",
    # Computation
    "compute_detailed_results",
    "compute_permuted_results",
    # Data
    "read_grp",
    "read_gmt",
    "load_rank_matrix",
    "load_cmap_data",
    # Types
    "CmapData",
    "DetailedResults",
    "PermutedResults",
]


def cmap(
    up_genes: list[str] | str | Path,
    down_genes: list[str] | str | Path,
    data_dir: str | Path = ".",
) -> tuple[DetailedResults, PermutedResults, PermutedResults]:
    """Run the full CMAP analysis pipeline.

    This is the high-level entry point mirroring ``cmap.default()`` from
    the original R implementation.  It computes:

    1. **Detailed results** — per-instance KS scores and connectivity.
    2. **Permuted results by name** — aggregated by drug name.
    3. **Permuted results by cell** — aggregated by drug name + cell line.

    Parameters
    ----------
    up_genes:
        Up-regulated gene identifiers.  Accepts a ``list[str]``, a path to
        a ``.grp`` file (str or Path).
    down_genes:
        Down-regulated gene identifiers (same input types as *up_genes*).
    data_dir:
        Directory containing the CMAP data files: rank matrix, ``.data``
        metadata files, and optionally pre-computed ``.npz`` files for
        p-values and specificity.

    Returns
    -------
    tuple[DetailedResults, PermutedResults, PermutedResults]
        ``(detailed_results, permuted_by_name, permuted_by_cell)``

    Examples
    --------
    >>> det, by_name, by_cell = cmap("readUp.grp", "readDown.grp", "data/")
    >>> det.to_dataframe().head(10)
    >>> by_name.to_dataframe().loc["geldanamycin"]
    >>> by_cell.to_dataframe().loc["tanespimycin - MCF7"]
    """
    data_dir = Path(data_dir)

    # Load all CMAP data
    data = load_cmap_data(data_dir)

    # Load pre-computed p-value distributions if available
    pv_path = data_dir / "pValue.npz"
    p_value_dist = load_p_values(pv_path) if pv_path.exists() else {}

    pvc_path = data_dir / "pValueByCell.npz"
    p_value_cell_dist = load_p_values(pvc_path) if pvc_path.exists() else {}

    # Load pre-computed specificity matrices if available
    sp_path = data_dir / "specificity.npz"
    spec_by_name = load_specificity(sp_path) if sp_path.exists() else None

    spc_path = data_dir / "specificityCell.npz"
    spec_by_cell = load_specificity(spc_path) if spc_path.exists() else None

    # Compute detailed results
    det = compute_detailed_results(up_genes, down_genes, data)

    # Compute permuted results by name
    by_name = compute_permuted_results(
        det, data, mode="by_name",
        p_value_dist=p_value_dist,
        specificity_matrix=spec_by_name,
    )

    # Compute permuted results by cell line
    by_cell = compute_permuted_results(
        det, data, mode="by_cell",
        p_value_dist=p_value_cell_dist,
        specificity_matrix=spec_by_cell,
    )

    return det, by_name, by_cell
