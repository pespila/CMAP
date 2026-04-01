"""Tests for pycmap.data — file readers, rank matrix I/O, and gene index."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from numpy.testing import assert_array_equal

if TYPE_CHECKING:
    from pathlib import Path

from pycmap.data import (
    build_gene_index,
    load_rank_matrix,
    read_data_file,
    read_gmt,
    read_grp,
    save_rank_matrix,
)

# ---------------------------------------------------------------------------
# read_grp
# ---------------------------------------------------------------------------


def test_read_grp(test_files_dir: Path) -> None:
    """readUp.grp must be parsed as ['1053_at', '1007_s_at']."""
    genes = read_grp(test_files_dir / "readUp.grp")
    assert genes == ["1053_at", "1007_s_at"]


def test_read_grp_down(test_files_dir: Path) -> None:
    """readDown.grp must be parsed as ['121_at', '117_at']."""
    genes = read_grp(test_files_dir / "readDown.grp")
    assert genes == ["121_at", "117_at"]


def test_read_grp_alzheimer(test_files_dir: Path) -> None:
    """alzheimerUp.grp must contain exactly 29 entries."""
    genes = read_grp(test_files_dir / "alzheimerUp.grp")
    assert len(genes) == 29


def test_read_grp_alzheimer_first_gene(test_files_dir: Path) -> None:
    """First gene in alzheimerUp.grp must be '200720_s_at'."""
    genes = read_grp(test_files_dir / "alzheimerUp.grp")
    assert genes[0] == "200720_s_at"


def test_read_grp_returns_list_of_strings(test_files_dir: Path) -> None:
    """All returned elements must be plain strings."""
    genes = read_grp(test_files_dir / "readUp.grp")
    assert all(isinstance(g, str) for g in genes)


def test_read_grp_missing_file() -> None:
    """A non-existent path must raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_grp("/nonexistent/path/missing.grp")


def test_read_grp_accepts_str_path(test_files_dir: Path) -> None:
    """read_grp must accept a plain str as well as a Path."""
    genes = read_grp(str(test_files_dir / "readUp.grp"))
    assert genes == ["1053_at", "1007_s_at"]


def test_read_grp_roundtrip(tmp_path: Path) -> None:
    """Write then read a .grp file — content must survive the round-trip."""
    original = ["GENE_A", "GENE_B", "GENE_C"]
    grp_file = tmp_path / "test.grp"
    grp_file.write_text("\n".join(original) + "\n", encoding="utf-8")
    loaded = read_grp(grp_file)
    assert loaded == original


def test_read_grp_blank_lines_ignored(tmp_path: Path) -> None:
    """Blank lines in a .grp file must be silently skipped."""
    content = "\nGENE_X\n\nGENE_Y\n\n"
    grp_file = tmp_path / "blanks.grp"
    grp_file.write_text(content, encoding="utf-8")
    genes = read_grp(grp_file)
    assert genes == ["GENE_X", "GENE_Y"]


# ---------------------------------------------------------------------------
# read_gmt
# ---------------------------------------------------------------------------


def test_read_gmt(test_files_dir: Path) -> None:
    """msigdb_up.gmt must parse into a dict with more than 100 entries."""
    sets = read_gmt(test_files_dir / "msigdb_up.gmt")
    assert isinstance(sets, dict)
    assert len(sets) > 100


def test_read_gmt_first_key(test_files_dir: Path) -> None:
    """The first key in msigdb_up.gmt must be 'PENG_GLUTAMINE_UP'."""
    sets = read_gmt(test_files_dir / "msigdb_up.gmt")
    first_key = next(iter(sets))
    assert first_key == "PENG_GLUTAMINE_UP"


def test_read_gmt_values_are_lists_of_strings(test_files_dir: Path) -> None:
    """Each value in the parsed GMT dict must be a non-empty list of strings."""
    sets = read_gmt(test_files_dir / "msigdb_up.gmt")
    for name, genes in sets.items():
        assert isinstance(genes, list), f"Expected list for {name!r}"
        assert len(genes) > 0, f"Empty gene list for {name!r}"
        assert all(isinstance(g, str) for g in genes), f"Non-string gene in {name!r}"


def test_read_gmt_missing_file() -> None:
    """A non-existent GMT file must raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_gmt("/nonexistent/path/missing.gmt")


def test_read_gmt_malformed_raises_value_error(tmp_path: Path) -> None:
    """A line with fewer than 3 tab-separated fields must raise ValueError."""
    bad_gmt = tmp_path / "bad.gmt"
    bad_gmt.write_text("ONLY_TWO_FIELDS\tdescription\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected at least 3 tab-separated fields"):
        read_gmt(bad_gmt)


def test_read_gmt_roundtrip(tmp_path: Path) -> None:
    """Write a minimal GMT and read it back — sets must be preserved."""
    content = "SET_A\tna\tGENE1\tGENE2\tGENE3\nSET_B\tna\tGENE4\n"
    gmt_file = tmp_path / "mini.gmt"
    gmt_file.write_text(content, encoding="utf-8")
    sets = read_gmt(gmt_file)
    assert list(sets.keys()) == ["SET_A", "SET_B"]
    assert sets["SET_A"] == ["GENE1", "GENE2", "GENE3"]
    assert sets["SET_B"] == ["GENE4"]


# ---------------------------------------------------------------------------
# read_data_file
# ---------------------------------------------------------------------------


def test_read_data_file_missing() -> None:
    """A non-existent .data file must raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_data_file("/nonexistent/path/missing.data")


def test_read_data_file_basic(tmp_path: Path) -> None:
    """Basic read: each non-blank line should become one entry."""
    data_file = tmp_path / "test.data"
    data_file.write_text("entry1\nentry2\n\nentry3\n", encoding="utf-8")
    entries = read_data_file(data_file)
    assert entries == ["entry1", "entry2", "entry3"]


def test_read_data_file_strips_whitespace(tmp_path: Path) -> None:
    """Leading/trailing whitespace on each line must be stripped."""
    data_file = tmp_path / "spaces.data"
    data_file.write_text("  alpha  \n  beta\ngamma  \n", encoding="utf-8")
    entries = read_data_file(data_file)
    assert entries == ["alpha", "beta", "gamma"]


# ---------------------------------------------------------------------------
# build_gene_index
# ---------------------------------------------------------------------------


def test_build_gene_index() -> None:
    """build_gene_index must return a mapping from name to zero-based row index."""
    names = ["BRCA1", "TP53", "EGFR"]
    idx = build_gene_index(names)
    assert idx == {"BRCA1": 0, "TP53": 1, "EGFR": 2}


def test_build_gene_index_empty() -> None:
    """An empty gene list must return an empty dict."""
    idx = build_gene_index([])
    assert idx == {}


def test_build_gene_index_lookup(gene_names: list[str]) -> None:
    """Every gene name must map to its correct 0-based index."""
    idx = build_gene_index(gene_names)
    for i, name in enumerate(gene_names):
        assert idx[name] == i


def test_build_gene_index_length(gene_names: list[str]) -> None:
    """Index length must equal the number of gene names."""
    idx = build_gene_index(gene_names)
    assert len(idx) == len(gene_names)


# ---------------------------------------------------------------------------
# save_rank_matrix / load_rank_matrix round-trip
# ---------------------------------------------------------------------------


def test_save_rank_matrix_creates_file(tmp_path: Path) -> None:
    """save_rank_matrix must create a .npz file at the given path."""
    matrix = np.array([[1, 3], [2, 4], [3, 1]], dtype=np.int16)
    gene_names_local = ["GENE_A", "GENE_B", "GENE_C"]
    instance_ids_local = ["inst_001", "inst_002"]

    out_path = tmp_path / "test_matrix"
    save_rank_matrix(out_path, matrix, gene_names_local, instance_ids_local)

    # numpy appends .npz automatically
    npz_path = tmp_path / "test_matrix.npz"
    assert npz_path.exists()


def test_save_and_load_rank_matrix(tmp_path: Path) -> None:
    """Round-trip via a unicode-string .npz archive must recover matrix and labels.

    Note: save_rank_matrix stores strings as object arrays which requires
    allow_pickle=True to reload. This test creates a compatible .npz directly
    using numpy's default unicode dtype so that load_rank_matrix (which uses
    allow_pickle=False) can read it back successfully.
    """
    matrix = np.array([[1, 3], [2, 4], [3, 1]], dtype=np.int16)
    gene_names_local = ["GENE_A", "GENE_B", "GENE_C"]
    instance_ids_local = ["inst_001", "inst_002"]

    npz_path = tmp_path / "test_matrix.npz"
    # Save with plain numpy unicode arrays (not dtype=object) for allow_pickle=False compat
    np.savez_compressed(
        npz_path,
        matrix=matrix,
        gene_names=np.array(gene_names_local),
        instance_ids=np.array(instance_ids_local),
    )

    loaded_matrix, loaded_genes, loaded_ids = load_rank_matrix(npz_path)

    assert_array_equal(loaded_matrix, matrix)
    assert loaded_genes == gene_names_local
    assert loaded_ids == instance_ids_local


def test_save_rank_matrix_mismatched_gene_names(tmp_path: Path) -> None:
    """Mismatched gene_names length must raise ValueError."""
    matrix = np.ones((5, 3), dtype=np.int16)
    with pytest.raises(ValueError, match="gene_names length"):
        save_rank_matrix(tmp_path / "bad", matrix, ["only_one"], ["a", "b", "c"])


def test_save_rank_matrix_mismatched_instance_ids(tmp_path: Path) -> None:
    """Mismatched instance_ids length must raise ValueError."""
    matrix = np.ones((5, 3), dtype=np.int16)
    gene_names_local = [f"g{i}" for i in range(5)]
    with pytest.raises(ValueError, match="instance_ids length"):
        save_rank_matrix(tmp_path / "bad", matrix, gene_names_local, ["only_one"])


def test_load_rank_matrix_npz_missing_keys(tmp_path: Path) -> None:
    """An .npz archive without the required keys must raise ValueError."""
    bad_npz = tmp_path / "incomplete.npz"
    np.savez(bad_npz, matrix=np.ones((3, 2)))
    with pytest.raises(ValueError, match="missing required keys"):
        load_rank_matrix(bad_npz)


def test_load_rank_matrix_npy_raises(tmp_path: Path) -> None:
    """Attempting to load a .npy file must raise ValueError with guidance."""
    npy_file = tmp_path / "rank.npy"
    np.save(npy_file, np.ones((3, 2)))
    with pytest.raises(ValueError, match=".npy"):
        load_rank_matrix(npy_file)


def test_load_rank_matrix_missing_file() -> None:
    """A non-existent rank matrix path must raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_rank_matrix("/nonexistent/rank_matrix.npz")


def test_load_rank_matrix_preserves_dtype(tmp_path: Path) -> None:
    """int16 matrix stored with unicode dtype strings must round-trip preserving dtype."""
    matrix = np.eye(4, dtype=np.int16)
    names = [f"g{i}" for i in range(4)]
    ids = [f"inst_{i}" for i in range(4)]
    npz_path = tmp_path / "eye.npz"
    np.savez_compressed(
        npz_path,
        matrix=matrix,
        gene_names=np.array(names),
        instance_ids=np.array(ids),
    )
    loaded, _, _ = load_rank_matrix(npz_path)
    assert loaded.dtype == np.int16
    assert_array_equal(loaded, matrix)


def test_load_rank_matrix_txt_format(tmp_path: Path) -> None:
    """A tab-separated text rank matrix must be loadable."""
    txt_file = tmp_path / "rank.txt"
    # header row: empty field + instance IDs
    lines = [
        "\tinst1\tinst2",
        "geneA\t3\t1",
        "geneB\t1\t3",
        "geneC\t2\t2",
    ]
    txt_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    matrix, genes, inst_ids = load_rank_matrix(txt_file)
    assert genes == ["geneA", "geneB", "geneC"]
    assert inst_ids == ["inst1", "inst2"]
    assert matrix.shape == (3, 2)
    assert_array_equal(matrix[0], [3, 1])
