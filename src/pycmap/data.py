"""Data loading for .grp, .gmt, rank matrix, and .data files.

This module provides all file I/O required by the CMAP pipeline.  It supports
the legacy plain-text formats used by the original R implementation as well as
NumPy binary formats (.npy / .npz) introduced by this Python port.

Design decisions
----------------
- No pandas — all parsing is done with plain Python and NumPy for minimal
  import cost and maximum throughput on large rank matrices.
- :class:`pathlib.Path` is used throughout for OS-agnostic path handling.
- All public functions raise :exc:`FileNotFoundError` or :exc:`ValueError` with
  descriptive messages so that callers get actionable diagnostics.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pycmap._types import CmapData

# ---------------------------------------------------------------------------
# Gene-set file readers
# ---------------------------------------------------------------------------


def read_grp(path: str | Path) -> list[str]:
    """Read a ``.grp`` file containing one gene identifier per line.

    Whitespace is stripped from every line; blank lines are silently skipped.

    Parameters
    ----------
    path:
        Path to the ``.grp`` file.

    Returns
    -------
    list[str]
        Ordered list of gene identifiers.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.

    Examples
    --------
    >>> genes = read_grp("readUp.grp")
    >>> genes[0]
    '1053_at'
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"GRP file not found: {p}")

    return [line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read a GMT (Gene Matrix Transposed) file.

    Each line is tab-separated with the layout::

        <set_name>  <description_or_url>  <gene1>  <gene2>  ...

    The description / URL field (column index 1) is ignored.

    Parameters
    ----------
    path:
        Path to the ``.gmt`` file.

    Returns
    -------
    dict[str, list[str]]
        Mapping from gene-set name to the list of gene identifiers in that set.
        Order of the sets reflects their order in the file.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    ValueError
        When a line has fewer than three tab-separated fields (i.e., the set
        has no genes).

    Examples
    --------
    >>> sets = read_gmt("msigdb_up.gmt")
    >>> list(sets.keys())[0]
    'PENG_GLUTAMINE_UP'
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"GMT file not found: {p}")

    result: dict[str, list[str]] = {}
    for lineno, raw in enumerate(p.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError(
                f"GMT file {p!r}, line {lineno}: expected at least 3 tab-separated fields, "
                f"got {len(fields)}"
            )
        set_name = fields[0]
        # fields[1] is description/URL — intentionally ignored
        genes = [g for g in fields[2:] if g]
        result[set_name] = genes

    return result


# ---------------------------------------------------------------------------
# Generic .data file reader
# ---------------------------------------------------------------------------


def read_data_file(path: str | Path) -> list[str]:
    """Read a ``.data`` file containing one entry per line.

    Used for instance IDs, instance names, cell lines, and similar per-column
    metadata.  Whitespace is stripped; blank lines are skipped.

    Parameters
    ----------
    path:
        Path to the ``.data`` file.

    Returns
    -------
    list[str]
        Ordered list of string entries.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Data file not found: {p}")

    return [line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# Rank matrix loaders
# ---------------------------------------------------------------------------


def _load_rank_matrix_npz(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    """Load a rank matrix from a NumPy ``.npz`` archive."""
    archive = np.load(path, allow_pickle=False)

    missing = {"matrix", "gene_names", "instance_ids"} - set(archive.files)
    if missing:
        raise ValueError(
            f"NPZ file {path!r} is missing required keys: {sorted(missing)}. "
            "Expected keys: 'matrix', 'gene_names', 'instance_ids'."
        )

    matrix: np.ndarray = archive["matrix"]
    gene_names: list[str] = archive["gene_names"].tolist()
    instance_ids: list[str] = archive["instance_ids"].tolist()

    return matrix, gene_names, instance_ids


def _load_rank_matrix_directory(directory: Path) -> tuple[np.ndarray, list[str], list[str]]:
    """Load a rank matrix from the legacy R split-text format.

    Expects two files inside *directory*:

    - ``rankMatrix.txt.headn1`` — tab-separated header row; first field is
      empty, remaining fields are instance IDs.
    - ``rankMatrix.txt.sed1d`` — tab-separated data rows; first field is the
      gene name, remaining fields are integer ranks.
    """
    head_path = directory / "rankMatrix.txt.headn1"
    data_path = directory / "rankMatrix.txt.sed1d"

    if not head_path.exists():
        raise FileNotFoundError(
            f"Legacy rank matrix header not found: {head_path}. "
            "Expected 'rankMatrix.txt.headn1' inside the data directory."
        )
    if not data_path.exists():
        raise FileNotFoundError(
            f"Legacy rank matrix data not found: {data_path}. "
            "Expected 'rankMatrix.txt.sed1d' inside the data directory."
        )

    # Parse header: first field is an empty placeholder, rest are instance IDs
    header_line = head_path.read_text(encoding="utf-8").rstrip("\n")
    header_fields = header_line.split("\t")
    instance_ids: list[str] = header_fields[1:]  # drop the leading empty field

    n_instances = len(instance_ids)

    gene_names: list[str] = []
    rows: list[np.ndarray] = []

    for lineno, raw in enumerate(data_path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip():
            continue
        fields = raw.split("\t")
        if len(fields) != n_instances + 1:
            raise ValueError(
                f"{data_path!r}, line {lineno}: expected {n_instances + 1} tab-separated "
                f"fields (1 gene name + {n_instances} ranks), got {len(fields)}."
            )
        gene_names.append(fields[0])
        rows.append(np.array(fields[1:], dtype=np.int32))

    if not rows:
        raise ValueError(f"Rank matrix data file is empty: {data_path}")

    matrix = np.stack(rows, axis=0)  # (n_genes, n_instances)

    # Downcast to int16 when values fit — saves ~50 % memory for typical CMAP data
    if matrix.max() <= np.iinfo(np.int16).max:
        matrix = matrix.astype(np.int16)

    return matrix, gene_names, instance_ids


def _load_rank_matrix_txt(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    """Load a rank matrix from a single tab-separated text file.

    The file is expected to have a header row (first field empty or a label,
    remaining fields are instance IDs) and data rows where the first field is
    the gene name and the remaining fields are integer ranks.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    non_empty = [ln for ln in lines if ln.strip()]

    if not non_empty:
        raise ValueError(f"Rank matrix text file is empty: {path}")

    header_fields = non_empty[0].split("\t")
    instance_ids: list[str] = header_fields[1:]
    n_instances = len(instance_ids)

    gene_names: list[str] = []
    rows: list[np.ndarray] = []

    for lineno, raw in enumerate(non_empty[1:], start=2):
        fields = raw.split("\t")
        if len(fields) != n_instances + 1:
            raise ValueError(
                f"{path!r}, line {lineno}: expected {n_instances + 1} tab-separated "
                f"fields (1 gene name + {n_instances} ranks), got {len(fields)}."
            )
        gene_names.append(fields[0])
        rows.append(np.array(fields[1:], dtype=np.int32))

    if not rows:
        raise ValueError(f"Rank matrix text file has a header but no data rows: {path}")

    matrix = np.stack(rows, axis=0)  # (n_genes, n_instances)

    if matrix.max() <= np.iinfo(np.int16).max:
        matrix = matrix.astype(np.int16)

    return matrix, gene_names, instance_ids


def load_rank_matrix(path: str | Path) -> tuple[np.ndarray, list[str], list[str]]:
    """Load a rank matrix from one of the supported formats.

    Dispatch is based on the nature of *path*:

    +---------------------+------------------------------------------------+
    | Condition           | Action                                         |
    +=====================+================================================+
    | ``.npz`` file       | Load NumPy archive (keys: matrix, gene_names,  |
    |                     | instance_ids)                                  |
    +---------------------+------------------------------------------------+
    | ``.npy`` file       | Raise :exc:`ValueError` — use ``.npz`` instead |
    +---------------------+------------------------------------------------+
    | Directory           | Legacy R split-text format (headn1 + sed1d)    |
    +---------------------+------------------------------------------------+
    | ``.txt`` (or other) | Single tab-separated text file                 |
    +---------------------+------------------------------------------------+

    Parameters
    ----------
    path:
        Path to the rank matrix source (file or directory).

    Returns
    -------
    tuple[np.ndarray, list[str], list[str]]
        A 3-tuple of:

        - ``matrix`` — integer ndarray of shape ``(n_genes, n_instances)``
          with dtype ``int16`` or ``int32``.
        - ``gene_names`` — row labels, length ``n_genes``.
        - ``instance_ids`` — column labels, length ``n_instances``.

    Raises
    ------
    FileNotFoundError
        When *path* (or expected sub-files) does not exist.
    ValueError
        When a ``.npy`` file is provided, or when the file content is
        structurally invalid.
    """
    p = Path(path)

    if p.is_dir():
        return _load_rank_matrix_directory(p)

    if not p.exists():
        raise FileNotFoundError(f"Rank matrix path not found: {p}")

    suffix = p.suffix.lower()

    if suffix == ".npy":
        raise ValueError(
            f"Plain .npy files do not store gene or instance names.  "
            f"Save the rank matrix as .npz using save_rank_matrix() and load that instead: {p}"
        )

    if suffix == ".npz":
        return _load_rank_matrix_npz(p)

    # Fall-through: treat as tab-separated text
    return _load_rank_matrix_txt(p)


# ---------------------------------------------------------------------------
# Gene index helper
# ---------------------------------------------------------------------------


def build_gene_index(gene_names: list[str]) -> dict[str, int]:
    """Build a name-to-row-index mapping for O(1) gene lookups.

    Parameters
    ----------
    gene_names:
        Ordered list of gene identifiers as returned by :func:`load_rank_matrix`.

    Returns
    -------
    dict[str, int]
        Mapping from gene name to its zero-based row index in the rank matrix.

    Examples
    --------
    >>> idx = build_gene_index(["BRCA1", "TP53", "EGFR"])
    >>> idx["TP53"]
    1
    """
    return {name: i for i, name in enumerate(gene_names)}


# ---------------------------------------------------------------------------
# Full CMAP data loader
# ---------------------------------------------------------------------------


def load_cmap_data(data_dir: str | Path) -> CmapData:
    """Load all CMAP data from a directory and return a :class:`~pycmap._types.CmapData`.

    Expected directory contents
    ---------------------------
    The rank matrix may be present in any format supported by
    :func:`load_rank_matrix`:

    - ``rankMatrix.npz`` (preferred — fastest load)
    - ``rankMatrix.npy`` (not supported — raises an error)
    - The directory itself as the legacy split-text format (headn1 + sed1d)
    - ``rankMatrix.txt`` (single tab-separated text file)

    Additionally the following ``.data`` files must be present:

    - ``instID.data`` — numeric instance IDs (one per column)
    - ``instanceList.data`` — instance names per column
    - ``cellLineList.data`` — cell lines per column
    - ``instanceNames.data`` — unique instance names
    - ``instanceCellNames.data`` — unique ``"name - cellline"`` combinations

    Parameters
    ----------
    data_dir:
        Directory containing the CMAP data files.

    Returns
    -------
    CmapData
        Fully populated data container ready for scoring.

    Raises
    ------
    FileNotFoundError
        When *data_dir* does not exist or a required ``.data`` file is missing.
    ValueError
        When the rank matrix cannot be loaded (see :func:`load_rank_matrix`).
    """
    d = Path(data_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"CMAP data directory not found: {d}")

    # ------------------------------------------------------------------
    # 1. Rank matrix — probe for supported file formats in preference order
    # ------------------------------------------------------------------
    matrix: np.ndarray
    gene_names: list[str]
    instance_id_strs: list[str]

    npz_path = d / "rankMatrix.npz"
    txt_path = d / "rankMatrix.txt"

    if npz_path.exists():
        matrix, gene_names, instance_id_strs = load_rank_matrix(npz_path)
    elif txt_path.exists():
        matrix, gene_names, instance_id_strs = load_rank_matrix(txt_path)
    else:
        # Try the directory itself as the legacy split-text root
        matrix, gene_names, instance_id_strs = load_rank_matrix(d)

    # ------------------------------------------------------------------
    # 2. Per-column metadata
    # ------------------------------------------------------------------
    instance_ids_raw = read_data_file(d / "instID.data")
    instance_ids = np.array(instance_ids_raw, dtype=np.int64)

    instance_names = read_data_file(d / "instanceList.data")
    cell_lines = read_data_file(d / "cellLineList.data")
    unique_names = read_data_file(d / "instanceNames.data")
    unique_name_cell = read_data_file(d / "instanceCellNames.data")

    return CmapData(
        rank_matrix=matrix,
        gene_names=gene_names,
        instance_ids=instance_ids,
        instance_names=instance_names,
        cell_lines=cell_lines,
        unique_names=unique_names,
        unique_name_cell=unique_name_cell,
    )


# ---------------------------------------------------------------------------
# Rank matrix serialiser
# ---------------------------------------------------------------------------


def save_rank_matrix(
    path: str | Path,
    matrix: np.ndarray,
    gene_names: list[str],
    instance_ids: list[str],
) -> None:
    """Save a rank matrix to a ``.npz`` archive for fast subsequent loading.

    The archive stores three arrays under fixed keys:

    - ``matrix`` — the rank matrix ndarray.
    - ``gene_names`` — 1-D array of gene name strings.
    - ``instance_ids`` — 1-D array of instance ID strings.

    Parameters
    ----------
    path:
        Destination path.  The ``.npz`` extension is appended automatically by
        :func:`numpy.savez_compressed` if not already present.
    matrix:
        Integer ndarray of shape ``(n_genes, n_instances)``.
    gene_names:
        Row labels, length must equal ``matrix.shape[0]``.
    instance_ids:
        Column labels, length must equal ``matrix.shape[1]``.

    Raises
    ------
    ValueError
        When the length of *gene_names* or *instance_ids* does not match the
        corresponding matrix dimension.

    Examples
    --------
    >>> import numpy as np
    >>> m = np.array([[1, 2], [3, 4]], dtype=np.int16)
    >>> save_rank_matrix("/tmp/test_matrix", m, ["GENE_A", "GENE_B"], ["inst1", "inst2"])
    """
    p = Path(path)
    n_genes, n_instances = matrix.shape

    if len(gene_names) != n_genes:
        raise ValueError(
            f"gene_names length ({len(gene_names)}) does not match matrix row count ({n_genes})."
        )
    if len(instance_ids) != n_instances:
        raise ValueError(
            f"instance_ids length ({len(instance_ids)}) does not match "
            f"matrix column count ({n_instances})."
        )

    np.savez_compressed(
        p,
        matrix=matrix,
        gene_names=np.array(gene_names, dtype=object),
        instance_ids=np.array(instance_ids, dtype=object),
    )
