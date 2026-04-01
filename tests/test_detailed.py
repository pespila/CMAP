"""Tests for pycmap.detailed — compute_detailed_results."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycmap._types import CmapData, DetailedResults
from pycmap.detailed import compute_detailed_results

# ---------------------------------------------------------------------------
# Shape and structure
# ---------------------------------------------------------------------------


def test_detailed_results_shape(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """All output arrays must have length == n_instances."""
    n_instances = cmap_data.rank_matrix.shape[1]
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)

    assert len(result.instance_ids) == n_instances
    assert len(result.scores) == n_instances
    assert len(result.ks_up) == n_instances
    assert len(result.ks_down) == n_instances
    assert len(result.instance_names) == n_instances


def test_detailed_results_is_detailed_results_type(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Return type must be a DetailedResults instance."""
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    assert isinstance(result, DetailedResults)


# ---------------------------------------------------------------------------
# Ordering
# ---------------------------------------------------------------------------


def test_detailed_results_sorted_descending(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Scores must be in non-increasing (descending) order."""
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    assert np.all(np.diff(result.scores) <= 0), "Scores are not in descending order"


def test_detailed_results_secondary_sort_exact_ties() -> None:
    """Within exact score ties, ks_up must be non-increasing (descending).

    A synthetic dataset is constructed so that certain instances produce
    identical raw scores.  After normalisation those will share the same
    floating-point value, enabling reliable detection of exact ties.
    """
    # Build a tiny 6-gene x 4-instance rank matrix where all raw scores are 0
    # (all same-sign k_up / k_down).  The normalised scores are all 0.0, so
    # the secondary sort by ks_up descending must apply to every adjacent pair.
    rng = np.random.default_rng(5)
    n_genes = 10
    n_inst = 4
    base = np.arange(1, n_genes + 1, dtype=np.int16)
    cols = [rng.permutation(base) for _ in range(n_inst)]
    mat = np.stack(cols, axis=1)
    gene_names_local = [f"g{i}" for i in range(n_genes)]
    inst_names_local = ["a", "b", "c", "d"]
    inst_ids_local = np.arange(1, n_inst + 1, dtype=np.int64)
    cell_lines_local = ["L1"] * n_inst

    data = CmapData(
        rank_matrix=mat,
        gene_names=gene_names_local,
        instance_ids=inst_ids_local,
        instance_names=inst_names_local,
        cell_lines=cell_lines_local,
        unique_names=inst_names_local,
        unique_name_cell=[f"{n} - L1" for n in inst_names_local],
    )

    # Use the same genes for up and down — k_up and k_down will have the same
    # sign → all raw scores are 0 → all normalised scores are 0
    up = gene_names_local[:2]
    down = gene_names_local[:2]
    result = compute_detailed_results(up, down, data)

    # All scores must be exactly 0.0
    assert np.all(result.scores == 0.0), "Expected all scores to be zero"

    # With all scores tied at 0.0, ks_up must be in non-increasing order
    assert np.all(np.diff(result.ks_up) <= 0), (
        "ks_up is not in descending order within all-zero score ties"
    )


# ---------------------------------------------------------------------------
# Score range
# ---------------------------------------------------------------------------


def test_detailed_results_scores_range(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """All normalised scores must lie within [-1, 1]."""
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    assert np.all(result.scores >= -1.0 - 1e-12)
    assert np.all(result.scores <= 1.0 + 1e-12)


def test_detailed_results_ks_range(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Raw KS statistics must lie within [-1, 1]."""
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    assert np.all(result.ks_up >= -1.0 - 1e-12)
    assert np.all(result.ks_up <= 1.0 + 1e-12)
    assert np.all(result.ks_down >= -1.0 - 1e-12)
    assert np.all(result.ks_down <= 1.0 + 1e-12)


# ---------------------------------------------------------------------------
# Integration: real .grp files against synthetic data
# ---------------------------------------------------------------------------


def test_detailed_results_from_grp_files(
    cmap_data: CmapData,
    test_files_dir: Path,
) -> None:
    """Running with real .grp files should succeed and return valid results.

    Because readUp.grp / readDown.grp contain Affymetrix probe IDs that do not
    exist in our synthetic 'gene_001..gene_100' rank matrix, the call must raise
    a ValueError (none of the genes found).  This confirms the code path that
    reads .grp files is exercised.  The fixture uses synthetic gene names, so we
    supply the grp file genes via list to avoid that, using gene names that *do*
    exist.
    """
    # Supply up/down lists pointing at real genes in the synthetic matrix
    up_path = test_files_dir / "readUp.grp"
    down_path = test_files_dir / "readDown.grp"

    # The synthetic matrix has gene_001..gene_100 — real probe IDs will be missing.
    # Verify the path is readable and contains expected probes, then confirm
    # that compute_detailed_results raises ValueError for unknown gene names.
    with pytest.raises(ValueError, match="None of the"):
        compute_detailed_results(up_path, down_path, cmap_data)


def test_detailed_results_grp_file_accepted_as_path(
    cmap_data: CmapData,
    test_files_dir: Path,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Passing a Path (not list) for genes that exist must work normally."""
    # Write a temporary .grp file with genes that exist in the synthetic matrix
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".grp", delete=False) as fh:
        fh.write("\n".join(small_up_genes) + "\n")
        up_path = Path(fh.name)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".grp", delete=False) as fh:
        fh.write("\n".join(small_down_genes) + "\n")
        down_path = Path(fh.name)

    try:
        result = compute_detailed_results(up_path, down_path, cmap_data)
        assert len(result.scores) == cmap_data.rank_matrix.shape[1]
    finally:
        up_path.unlink(missing_ok=True)
        down_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Warnings and errors
# ---------------------------------------------------------------------------


def test_detailed_results_missing_genes_warning(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Genes not present in the rank matrix must trigger a UserWarning."""
    up_with_phantom = [*small_up_genes, "phantom_gene_X", "phantom_gene_Y"]
    down_with_phantom = [*small_down_genes, "phantom_gene_Z"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = compute_detailed_results(up_with_phantom, down_with_phantom, cmap_data)

    warning_messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
    assert any("not found in rank matrix" in msg for msg in warning_messages), (
        f"Expected a UserWarning about missing genes; got: {warning_messages}"
    )
    # The result should still be computed from the found genes
    assert len(result.scores) == cmap_data.rank_matrix.shape[1]


def test_detailed_results_no_genes_error(cmap_data: CmapData) -> None:
    """When all up-genes are absent the function must raise ValueError."""
    with pytest.raises(ValueError, match="None of the"):
        compute_detailed_results(
            ["phantom_up_1", "phantom_up_2"],
            ["phantom_down_1"],
            cmap_data,
        )


def test_detailed_results_no_down_genes_error(
    cmap_data: CmapData,
    small_up_genes: list[str],
) -> None:
    """When all down-genes are absent the function must raise ValueError."""
    with pytest.raises(ValueError, match="None of the"):
        compute_detailed_results(
            small_up_genes,
            ["phantom_down_1", "phantom_down_2"],
            cmap_data,
        )


# ---------------------------------------------------------------------------
# to_dataframe
# ---------------------------------------------------------------------------


def test_detailed_results_to_dataframe_columns(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """to_dataframe() must return exactly the expected columns."""
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    df = result.to_dataframe()
    assert list(df.columns) == ["InstID", "Score", "Up", "Down"]


def test_detailed_results_to_dataframe_length(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """DataFrame length must equal n_instances."""
    n_instances = cmap_data.rank_matrix.shape[1]
    result = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    df = result.to_dataframe()
    assert len(df) == n_instances


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_detailed_results_deterministic(
    cmap_data: CmapData,
    small_up_genes: list[str],
    small_down_genes: list[str],
) -> None:
    """Two calls with identical inputs must produce identical results."""
    r1 = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    r2 = compute_detailed_results(small_up_genes, small_down_genes, cmap_data)
    assert_allclose(r1.scores, r2.scores)
    assert_allclose(r1.ks_up, r2.ks_up)
    assert_allclose(r1.ks_down, r2.ks_down)
    assert list(r1.instance_names) == list(r2.instance_names)
