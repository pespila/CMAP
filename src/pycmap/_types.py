"""Data container types for CMAP results.

All containers are plain :func:`dataclasses.dataclass` instances so that they
carry no hidden logic, remain trivially serialisable, and are easy to construct
in tests.  Heavy dependencies (pandas) are imported lazily under
``TYPE_CHECKING`` and inside the ``to_dataframe`` methods to keep the import
cost of the package low.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np  # noqa: TCH002

if TYPE_CHECKING:
    import pandas as pd


# ---------------------------------------------------------------------------
# Input data container
# ---------------------------------------------------------------------------


@dataclass
class CmapData:
    """All data required to run a CMAP analysis.

    This container is populated by the data-loading layer and passed unchanged
    into the scoring functions.

    Attributes
    ----------
    rank_matrix:
        Integer rank matrix of shape ``(n_genes, n_instances)``.  Values are
        1-indexed ranks (dtype ``int16`` or ``int32``).  Rows correspond to
        genes, columns to instances (drug treatments).
    gene_names:
        Row labels of ``rank_matrix``.  Length must equal ``n_genes``.
    instance_ids:
        1-D integer array of instance IDs.  Length must equal ``n_instances``.
    instance_names:
        Per-instance drug names aligned with ``rank_matrix`` columns.  Length
        must equal ``n_instances``.
    cell_lines:
        Per-instance cell-line labels aligned with ``rank_matrix`` columns.
        Length must equal ``n_instances``.
    unique_names:
        Deduplicated drug names preserving first-occurrence order (equivalent
        to ``instanceNames.data`` from the R installation).
    unique_name_cell:
        Deduplicated ``"name - cellline"`` combination strings preserving
        first-occurrence order (equivalent to ``instanceCellNames.data``).
    """

    rank_matrix: np.ndarray
    gene_names: list[str]
    instance_ids: np.ndarray
    instance_names: list[str]
    cell_lines: list[str]
    unique_names: list[str]
    unique_name_cell: list[str]


# ---------------------------------------------------------------------------
# Detailed results
# ---------------------------------------------------------------------------


@dataclass
class DetailedResults:
    """Output of the detailed CMAP results computation.

    Corresponds to the ``detRes`` matrix produced by ``cmap.detResults`` in
    ``legacy/cmap.R``.  Scores are normalised to ``[-1, 1]``.

    Attributes
    ----------
    instance_ids:
        1-D integer array of instance IDs.
    scores:
        1-D float array of normalised connectivity scores in ``[-1, 1]``.
    ks_up:
        1-D float array of signed KS statistics for the up-regulated gene set.
    ks_down:
        1-D float array of signed KS statistics for the down-regulated gene set.
    instance_names:
        Row labels for the result (drug name or ``"name - cellline"``).  Used
        as the DataFrame index.
    """

    instance_ids: np.ndarray
    scores: np.ndarray
    ks_up: np.ndarray
    ks_down: np.ndarray
    instance_names: list[str]

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to a :class:`pandas.DataFrame` with standard CMAP columns.

        Returns
        -------
        pandas.DataFrame
            Columns: ``InstID``, ``Score``, ``Up``, ``Down``.
            Index: ``instance_names``.

        Examples
        --------
        >>> import numpy as np
        >>> dr = DetailedResults(
        ...     instance_ids=np.array([1, 2]),
        ...     scores=np.array([0.9, -0.5]),
        ...     ks_up=np.array([0.4, -0.3]),
        ...     ks_down=np.array([-0.3, 0.2]),
        ...     instance_names=["drugA", "drugB"],
        ... )
        >>> df = dr.to_dataframe()
        >>> list(df.columns)
        ['InstID', 'Score', 'Up', 'Down']
        """
        import pandas as pd  # noqa: PLC0415

        return pd.DataFrame(
            {
                "InstID": self.instance_ids,
                "Score": self.scores,
                "Up": self.ks_up,
                "Down": self.ks_down,
            },
            index=self.instance_names,
        )


# ---------------------------------------------------------------------------
# Permuted results
# ---------------------------------------------------------------------------


@dataclass
class PermutedResults:
    """Output of the permuted CMAP results computation.

    Corresponds to the ``nameRes`` / ``cellRes`` matrices produced by
    ``cmap.permResults`` in ``legacy/cmap.R``.

    Attributes
    ----------
    mean:
        1-D float array of mean normalised scores across all instances of each
        drug name (or name-cell combination).
    n:
        1-D integer array of instance frequencies (how many rank-matrix columns
        belong to each name).
    enrichment:
        1-D float array of enrichment scores (KS statistic over the rank
        positions of the name's instances within the detailed results).
    p_value:
        1-D float array of empirical p-values derived from the permutation null
        distribution.
    specificity:
        1-D float array of specificity scores in ``[0, 1]``.
    non_null_pct:
        1-D float array of non-null percentages (see
        :func:`~pycmap.core.non_null_pct`).
    names:
        Row labels — unique drug names or ``"name - cellline"`` strings.  Used
        as the DataFrame index.
    """

    mean: np.ndarray
    n: np.ndarray
    enrichment: np.ndarray
    p_value: np.ndarray
    specificity: np.ndarray
    non_null_pct: np.ndarray
    names: list[str]

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to a :class:`pandas.DataFrame` with standard CMAP columns.

        Returns
        -------
        pandas.DataFrame
            Columns: ``mean``, ``n``, ``enrichment``, ``p``, ``specificity``,
            ``% non-null``.
            Index: ``names``.

        Examples
        --------
        >>> import numpy as np
        >>> pr = PermutedResults(
        ...     mean=np.array([0.3]),
        ...     n=np.array([5]),
        ...     enrichment=np.array([0.25]),
        ...     p_value=np.array([0.04]),
        ...     specificity=np.array([0.9]),
        ...     non_null_pct=np.array([80.0]),
        ...     names=["drugA"],
        ... )
        >>> df = pr.to_dataframe()
        >>> list(df.columns)
        ['mean', 'n', 'enrichment', 'p', 'specificity', '% non-null']
        """
        import pandas as pd  # noqa: PLC0415

        return pd.DataFrame(
            {
                "mean": self.mean,
                "n": self.n,
                "enrichment": self.enrichment,
                "p": self.p_value,
                "specificity": self.specificity,
                "% non-null": self.non_null_pct,
            },
            index=self.names,
        )
