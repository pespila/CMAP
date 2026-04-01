"""Pre-computation of p-value distributions and specificity matrices.

This module implements the equivalent of ``cmap.install`` from
``legacy/cmap.R``.  It handles:

1. Converting the legacy rank matrix text files to NumPy ``.npz`` format.
2. Generating ``.data`` metadata files from the CMAP instance inventory.
3. Pre-computing p-value null distributions (100 000 permutations per
   unique instance frequency).
4. Pre-computing the specificity matrix (enrichment of every drug against
   every MSigDB gene set).

All expensive operations are vectorised with NumPy.  The specificity
computation can be parallelised across CPU cores via
:class:`concurrent.futures.ProcessPoolExecutor`.
"""

from __future__ import annotations

import concurrent.futures
import logging
from pathlib import Path

import numpy as np

from pycmap.core import ks_permutation_batch
from pycmap.data import (
    load_cmap_data,
    read_gmt,
    read_grp,
    save_rank_matrix,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Rank matrix conversion
# ---------------------------------------------------------------------------


def build_rank_matrix(
    head_path: str | Path,
    data_path: str | Path,
    output_path: str | Path,
) -> None:
    """Convert legacy split-text rank matrix files to ``.npz`` format.

    Parameters
    ----------
    head_path:
        Path to ``rankMatrix.txt.headn1`` (tab-separated header row).
    data_path:
        Path to ``rankMatrix.txt.sed1d`` (tab-separated data rows).
    output_path:
        Destination ``.npz`` file path.
    """
    head_path = Path(head_path)
    data_path = Path(data_path)

    header_line = head_path.read_text(encoding="utf-8").rstrip("\n")
    header_fields = header_line.split("\t")
    instance_ids: list[str] = header_fields[1:]

    gene_names: list[str] = []
    rows: list[np.ndarray] = []

    for raw in data_path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        fields = raw.split("\t")
        gene_names.append(fields[0])
        rows.append(np.array(fields[1:], dtype=np.int32))

    matrix = np.stack(rows, axis=0)
    if matrix.max() <= np.iinfo(np.int16).max:
        matrix = matrix.astype(np.int16)

    save_rank_matrix(output_path, matrix, gene_names, instance_ids)
    logger.info("Rank matrix saved to %s (%s)", output_path, matrix.shape)


# ---------------------------------------------------------------------------
# Data file generation from CMAP instance inventory
# ---------------------------------------------------------------------------


def create_data_files(
    cmap_data_path: str | Path,
    n_instances: int,
    output_dir: str | Path,
) -> None:
    """Generate ``.data`` metadata files from the CMAP instance inventory.

    This mirrors the R ``createDataFiles`` function.  The inventory file
    (``cmapData.data``) is a tab-separated dump of the CMAP Excel file.

    Parameters
    ----------
    cmap_data_path:
        Path to ``cmapData.data`` (tab-separated instance inventory).
    n_instances:
        Number of instances (columns) in the rank matrix.  Only the first
        *n_instances* rows of the inventory are used.
    output_dir:
        Directory where ``.data`` files will be written.
    """
    cmap_data_path = Path(cmap_data_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lines = cmap_data_path.read_text(encoding="utf-8").splitlines()
    # Skip header row (same as R: splittedLines <- splittedLines[-1])
    data_lines = lines[1:]

    inst_ids: list[str] = []
    instance_names: list[str] = []
    cell_lines: list[str] = []

    for line in data_lines[:n_instances]:
        fields = line.split("\t")
        inst_ids.append(fields[0])
        instance_names.append(fields[2])
        cell_lines.append(fields[6])

    # Write per-instance files
    _write_lines(output_dir / "instID.data", inst_ids)
    _write_lines(output_dir / "instanceList.data", instance_names)
    _write_lines(output_dir / "cellLineList.data", cell_lines)

    # Unique instance names (preserving first-occurrence order)
    unique_names = list(dict.fromkeys(instance_names))
    _write_lines(output_dir / "instanceNames.data", unique_names)

    # Unique "name - cellline" combinations
    name_cell = [f"{n} - {c}" for n, c in zip(instance_names, cell_lines, strict=True)]
    unique_name_cell = list(dict.fromkeys(name_cell))
    _write_lines(output_dir / "instanceCellNames.data", unique_name_cell)
    _write_lines(output_dir / "instanceAndCellLineList.data", name_cell)

    logger.info(
        "Data files created in %s: %d instances, %d unique names, %d unique name+cell",
        output_dir,
        n_instances,
        len(unique_names),
        len(unique_name_cell),
    )


def _write_lines(path: Path, lines: list[str]) -> None:
    """Write a list of strings to a file, one per line."""
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# P-value null distributions
# ---------------------------------------------------------------------------


def precompute_p_values(
    n_instances: int,
    unique_frequencies: list[int],
    n_permutations: int = 100_000,
    seed: int | None = None,
    output_path: str | Path | None = None,
) -> dict[str, np.ndarray]:
    """Pre-compute null KS distributions for each unique instance frequency.

    For each frequency ``f`` in *unique_frequencies*, this generates
    *n_permutations* random KS scores by sampling ``f`` items from
    ``1..n_instances`` and computing the absolute KS statistic.

    Parameters
    ----------
    n_instances:
        Total number of instances (rank matrix columns).
    unique_frequencies:
        Sorted list of unique instance frequencies (how many times each
        drug appears).  Frequencies of 1 are excluded (sentinel p-value).
    n_permutations:
        Number of random permutations per frequency.  Default: 100 000.
    seed:
        Random seed for reproducibility.
    output_path:
        If given, save the result as a ``.npz`` file.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping from frequency string (e.g. ``"2"``, ``"15"``) to a 1-D
        array of absolute KS scores of length *n_permutations*.
    """
    rng = np.random.default_rng(seed)
    result: dict[str, np.ndarray] = {}

    # Filter out frequency 1 (always gets sentinel p-value)
    freqs = sorted(f for f in unique_frequencies if f > 1)

    for freq in freqs:
        logger.info("Computing null distribution for frequency %d ...", freq)
        null_dist = ks_permutation_batch(
            n_total=n_instances,
            sample_size=freq,
            n_permutations=n_permutations,
            rng=rng,
        )
        result[str(freq)] = null_dist

    if output_path is not None:
        np.savez_compressed(str(output_path), **result)
        logger.info("P-value distributions saved to %s", output_path)

    return result


def load_p_values(path: str | Path) -> dict[str, np.ndarray]:
    """Load pre-computed p-value distributions from a ``.npz`` file."""
    archive = np.load(str(path), allow_pickle=False)
    return {key: archive[key] for key in archive.files}


# ---------------------------------------------------------------------------
# Specificity matrix
# ---------------------------------------------------------------------------


def _compute_specificity_column(
    args: tuple[int, str, str, str],
) -> tuple[int, np.ndarray]:
    """Worker function for parallel specificity computation.

    Computes enrichment scores for one MSigDB gene set against all drugs.
    """
    col_idx, up_grp_path, down_grp_path, data_dir = args

    # Import here to avoid circular imports in worker processes
    from pycmap.detailed import compute_detailed_results  # noqa: PLC0415
    from pycmap.permuted import compute_permuted_results  # noqa: PLC0415

    data = load_cmap_data(data_dir)

    up_genes = read_grp(up_grp_path)
    down_genes = read_grp(down_grp_path)

    det = compute_detailed_results(up_genes, down_genes, data)
    perm = compute_permuted_results(det, data, mode="by_name")

    return col_idx, perm.enrichment


def precompute_specificity(
    data_dir: str | Path,
    msig_up_path: str | Path,
    msig_down_path: str | Path,
    mode: str = "by_name",
    n_workers: int = 1,
    output_path: str | Path | None = None,
) -> np.ndarray:
    """Compute the specificity matrix: enrichment of each drug vs each MSigDB set.

    This runs the full detailed+permuted pipeline for each of the ~312
    MSigDB gene sets.  Each gene set is independent, so the computation can
    be parallelised across CPU cores.

    Parameters
    ----------
    data_dir:
        Directory containing rank matrix and ``.data`` files.
    msig_up_path:
        Path to the MSigDB up-regulation ``.gmt`` file.
    msig_down_path:
        Path to the MSigDB down-regulation ``.gmt`` file.
    mode:
        ``"by_name"`` or ``"by_cell"``.
    n_workers:
        Number of parallel workers.  ``1`` for sequential (default).
    output_path:
        If given, save the result as a ``.npz`` file.

    Returns
    -------
    np.ndarray
        Specificity matrix of shape ``(n_unique_drugs, n_msigdb_sets)``.
    """
    data_dir = Path(data_dir)
    msig_up = read_gmt(msig_up_path)
    msig_down = read_gmt(msig_down_path)
    n_sets = len(msig_up)

    # Determine number of unique drugs for matrix shape
    data = load_cmap_data(data_dir)
    n_drugs = len(data.unique_names) if mode == "by_name" else len(data.unique_name_cell)

    # Write temporary .grp files for each MSigDB set
    tmp_dir = data_dir / "_msigdb_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    tasks: list[tuple[int, str, str, str]] = []
    up_names = list(msig_up.keys())
    down_names = list(msig_down.keys())

    for i in range(n_sets):
        up_file = tmp_dir / f"msigdb_up_{i}.grp"
        down_file = tmp_dir / f"msigdb_down_{i}.grp"
        _write_lines(up_file, msig_up[up_names[i]])
        _write_lines(down_file, msig_down[down_names[i]])
        tasks.append((i, str(up_file), str(down_file), str(data_dir)))

    spec_matrix = np.zeros((n_drugs, n_sets), dtype=np.float64)

    if n_workers <= 1:
        for task in tasks:
            col_idx, enrichments = _compute_specificity_column(task)
            spec_matrix[:, col_idx] = enrichments
            logger.info("Specificity: completed gene set %d/%d", col_idx + 1, n_sets)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_compute_specificity_column, task): task[0]
                for task in tasks
            }
            for future in concurrent.futures.as_completed(futures):
                col_idx, enrichments = future.result()
                spec_matrix[:, col_idx] = enrichments
                logger.info("Specificity: completed gene set %d/%d", col_idx + 1, n_sets)

    # Clean up temp files
    for f in tmp_dir.iterdir():
        f.unlink()
    tmp_dir.rmdir()

    if output_path is not None:
        np.savez_compressed(str(output_path), specificity=spec_matrix)
        logger.info("Specificity matrix saved to %s", output_path)

    return spec_matrix


def load_specificity(path: str | Path) -> np.ndarray:
    """Load a pre-computed specificity matrix from a ``.npz`` file."""
    archive = np.load(str(path), allow_pickle=False)
    return archive["specificity"]


# ---------------------------------------------------------------------------
# Full installation orchestrator
# ---------------------------------------------------------------------------


def install(
    data_dir: str | Path,
    n_permutations: int = 100_000,
    n_workers: int = 1,
    seed: int | None = None,
) -> None:
    """Run the full CMAP installation: pre-compute all required data files.

    This is the equivalent of ``cmap.install`` from ``legacy/cmap.R``.

    Prerequisites: the *data_dir* must already contain:
    - A rank matrix (in any format supported by :func:`~pycmap.data.load_rank_matrix`)
    - ``.data`` metadata files (or ``cmapData.data`` to generate them)
    - ``msigdb_up.gmt`` and ``msigdb_down.gmt``

    Parameters
    ----------
    data_dir:
        Directory containing all input data.
    n_permutations:
        Number of permutations for p-value null distributions.
    n_workers:
        Number of parallel workers for specificity computation.
    seed:
        Random seed for reproducibility.
    """
    data_dir = Path(data_dir)
    logger.info("Starting CMAP installation in %s", data_dir)

    # Load data to determine instance counts and frequencies
    data = load_cmap_data(data_dir)
    n_instances = data.rank_matrix.shape[1]

    # Compute unique frequencies for p-value distributions
    from collections import Counter  # noqa: PLC0415

    name_counts = Counter(data.instance_names)
    unique_freqs_by_name = sorted(set(name_counts.values()))

    name_cell = [
        f"{n} - {c}" for n, c in zip(data.instance_names, data.cell_lines, strict=True)
    ]
    cell_counts = Counter(name_cell)
    unique_freqs_by_cell = sorted(set(cell_counts.values()))

    # Pre-compute p-value distributions
    logger.info("Pre-computing p-value distributions (by name)...")
    precompute_p_values(
        n_instances, unique_freqs_by_name, n_permutations, seed,
        output_path=data_dir / "pValue.npz",
    )

    logger.info("Pre-computing p-value distributions (by cell)...")
    precompute_p_values(
        n_instances, unique_freqs_by_cell, n_permutations, seed,
        output_path=data_dir / "pValueByCell.npz",
    )

    # Pre-compute specificity matrices
    msig_up = data_dir / "msigdb_up.gmt"
    msig_down = data_dir / "msigdb_down.gmt"

    if msig_up.exists() and msig_down.exists():
        logger.info("Pre-computing specificity matrix (by name)...")
        precompute_specificity(
            data_dir, msig_up, msig_down, mode="by_name",
            n_workers=n_workers, output_path=data_dir / "specificity.npz",
        )

        logger.info("Pre-computing specificity matrix (by cell)...")
        precompute_specificity(
            data_dir, msig_up, msig_down, mode="by_cell",
            n_workers=n_workers, output_path=data_dir / "specificityCell.npz",
        )
    else:
        logger.warning(
            "MSigDB files not found (%s, %s). Skipping specificity pre-computation.",
            msig_up, msig_down,
        )

    logger.info("CMAP installation complete.")
