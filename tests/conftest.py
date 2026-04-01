"""Shared test fixtures for the pycmap-score test suite."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycmap._types import CmapData

# ---------------------------------------------------------------------------
# Constants used across fixtures
# ---------------------------------------------------------------------------

_N_GENES = 100
_N_INSTANCES = 20
_DRUG_NAMES = ["drugA", "drugB", "drugC", "drugD", "drugE"]
_CELL_LINES_POOL = ["MCF7", "PC3", "VCAP"]


# ---------------------------------------------------------------------------
# Rank-matrix fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def gene_names() -> list[str]:
    """100 gene name strings: gene_001 .. gene_100."""
    return [f"gene_{i:03d}" for i in range(1, _N_GENES + 1)]


@pytest.fixture(scope="session")
def instance_names() -> list[str]:
    """20 instance names — 5 unique drug names each repeated 4 times."""
    return [name for name in _DRUG_NAMES for _ in range(4)]


@pytest.fixture(scope="session")
def instance_ids() -> np.ndarray:
    """20 integer instance IDs (1-indexed)."""
    return np.arange(1, _N_INSTANCES + 1, dtype=np.int64)


@pytest.fixture(scope="session")
def cell_lines() -> list[str]:
    """20 cell-line strings drawn round-robin from 3 unique cell lines."""
    return [_CELL_LINES_POOL[i % len(_CELL_LINES_POOL)] for i in range(_N_INSTANCES)]


@pytest.fixture(scope="session")
def small_rank_matrix() -> np.ndarray:
    """100-gene x 20-instance int16 array where each column is a permutation of 1..100."""
    rng = np.random.default_rng(42)
    base = np.arange(1, _N_GENES + 1, dtype=np.int16)
    cols = [rng.permutation(base) for _ in range(_N_INSTANCES)]
    return np.stack(cols, axis=1)  # shape (100, 20)


@pytest.fixture(scope="session")
def cmap_data(
    small_rank_matrix: np.ndarray,
    gene_names: list[str],
    instance_names: list[str],
    instance_ids: np.ndarray,
    cell_lines: list[str],
) -> CmapData:
    """A fully assembled CmapData fixture built from the small synthetic matrix."""
    unique_names = list(dict.fromkeys(instance_names))  # preserve first-occurrence order
    unique_name_cell = list(
        dict.fromkeys(
            f"{name} - {cell}" for name, cell in zip(instance_names, cell_lines, strict=True)
        )
    )

    return CmapData(
        rank_matrix=small_rank_matrix,
        gene_names=gene_names,
        instance_ids=instance_ids,
        instance_names=instance_names,
        cell_lines=cell_lines,
        unique_names=unique_names,
        unique_name_cell=unique_name_cell,
    )


# ---------------------------------------------------------------------------
# Gene-set fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def small_up_genes(gene_names: list[str]) -> list[str]:
    """5 gene names that exist in the rank matrix (first 5)."""
    return gene_names[:5]


@pytest.fixture(scope="session")
def small_down_genes(gene_names: list[str]) -> list[str]:
    """5 gene names that exist in the rank matrix (last 5)."""
    return gene_names[-5:]


# ---------------------------------------------------------------------------
# Path fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def test_files_dir() -> Path:
    """Path to the test_files directory."""
    return Path("/Users/michael.wallner/Developer/CMAP/test_files")
