"""Data fetching from clue.io and legacy Broad Institute sources.

This module provides functions to download CMAP data from:

1. **clue.io** — the modern Connectivity Map platform (requires API key,
   register at https://clue.io).
2. **Legacy sources** — the original Broad Institute FTP archives and
   longevityTools data (no authentication required, but data is from the
   2006 release).

All network I/O uses ``httpx`` which is an **optional dependency**.
Install it with::

    pip install pycmap-score[fetch]

Functions
---------
- :func:`fetch_from_clue` — download data via the clue.io REST API.
- :func:`fetch_legacy_data` — download legacy rank matrix and metadata.
- :func:`fetch_msigdb_gene_sets` — download MSigDB .gmt files.
- :func:`fetch_and_install` — high-level: fetch data then run installation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)

_CLUE_API_BASE = "https://api.clue.io/api"
_LEGACY_RANK_MATRIX_URL = "https://portals.broadinstitute.org/cmap/rankMatrix.txt.zip"
_LONGEVITY_TOOLS_BASE = "https://girke.bioinformatics.ucr.edu/longevityTools"


def _ensure_httpx():  # type: ignore[no-untyped-def]
    """Import and return httpx, raising a helpful error if not installed."""
    try:
        import httpx  # noqa: PLC0415

        return httpx
    except ImportError as exc:
        raise ImportError(
            "httpx is required for data fetching. Install it with: pip install pycmap-score[fetch]"
        ) from exc


# ---------------------------------------------------------------------------
# clue.io API
# ---------------------------------------------------------------------------


def fetch_from_clue(
    api_key: str,
    output_dir: str | Path,
    *,
    base_url: str = _CLUE_API_BASE,
) -> Path:
    """Download CMAP data from the clue.io REST API.

    Requires registration at https://clue.io to obtain an API key.

    Parameters
    ----------
    api_key:
        Your clue.io API key (``user_key`` parameter).
    output_dir:
        Directory where downloaded files will be saved.
    base_url:
        Base URL for the clue.io API.  Override for testing.

    Returns
    -------
    Path
        The *output_dir* path containing the downloaded files.

    Raises
    ------
    ImportError
        When ``httpx`` is not installed.
    httpx.HTTPStatusError
        When the API returns an error response.
    """
    httpx = _ensure_httpx()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    headers = {"user_key": api_key, "Accept": "application/json"}

    with httpx.Client(base_url=base_url, headers=headers, timeout=120.0) as client:
        # Fetch instance metadata
        logger.info("Fetching instance metadata from clue.io...")
        resp = client.get("/perts", params={"q": '{"pert_type":"trt_cp"}', "l": 100})
        resp.raise_for_status()

        meta_path = output_dir / "clue_instances.json"
        meta_path.write_bytes(resp.content)
        logger.info("Instance metadata saved to %s (%d bytes)", meta_path, len(resp.content))

    logger.info("clue.io data download complete. Files in: %s", output_dir)
    return output_dir


# ---------------------------------------------------------------------------
# Legacy data sources
# ---------------------------------------------------------------------------


def fetch_legacy_data(
    output_dir: str | Path,
    *,
    rank_matrix_url: str = _LEGACY_RANK_MATRIX_URL,
) -> Path:
    """Download the legacy CMAP rank matrix from the Broad Institute.

    The original rank matrix (~500 MB compressed) is downloaded and
    extracted.  This data is from the 2006 CMAP release.

    Parameters
    ----------
    output_dir:
        Directory where downloaded files will be saved.
    rank_matrix_url:
        URL for the rank matrix ZIP file.

    Returns
    -------
    Path
        The *output_dir* path containing the extracted files.

    Raises
    ------
    ImportError
        When ``httpx`` is not installed.
    """
    import zipfile  # noqa: PLC0415

    httpx = _ensure_httpx()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    zip_path = output_dir / "rankMatrix.txt.zip"

    logger.info("Downloading legacy rank matrix from %s ...", rank_matrix_url)
    with httpx.stream("GET", rank_matrix_url, timeout=600.0, follow_redirects=True) as resp:
        resp.raise_for_status()
        with zip_path.open("wb") as f:
            for chunk in resp.iter_bytes(chunk_size=8192):
                f.write(chunk)

    logger.info("Extracting %s ...", zip_path)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(output_dir)

    zip_path.unlink()
    logger.info("Legacy rank matrix extracted to %s", output_dir)
    return output_dir


def fetch_msigdb_gene_sets(
    output_dir: str | Path,
    *,
    base_url: str = _LONGEVITY_TOOLS_BASE,
) -> tuple[Path, Path]:
    """Download MSigDB up/down gene set files.

    These ``.gmt`` files define the gene sets used for specificity
    computation in the CMAP algorithm.

    Parameters
    ----------
    output_dir:
        Directory where the ``.gmt`` files will be saved.
    base_url:
        Base URL for the longevityTools data.

    Returns
    -------
    tuple[Path, Path]
        Paths to ``msigdb_up.gmt`` and ``msigdb_down.gmt``.

    Raises
    ------
    ImportError
        When ``httpx`` is not installed.
    """
    httpx = _ensure_httpx()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []
    for name in ("msigdb_up.gmt", "msigdb_down.gmt"):
        url = f"{base_url}/data/{name}"
        logger.info("Downloading %s ...", url)
        resp = httpx.get(url, timeout=120.0, follow_redirects=True)
        resp.raise_for_status()
        dest = output_dir / name
        dest.write_bytes(resp.content)
        logger.info("Saved %s (%d bytes)", dest, len(resp.content))
        paths.append(dest)

    return paths[0], paths[1]


# ---------------------------------------------------------------------------
# High-level: fetch + install
# ---------------------------------------------------------------------------


def fetch_and_install(
    source: Literal["clue", "legacy", "both"] = "both",
    output_dir: str | Path = ".",
    *,
    api_key: str | None = None,
    n_permutations: int = 100_000,
    n_workers: int = 1,
    seed: int | None = None,
) -> Path:
    """Fetch CMAP data and run the full installation pipeline.

    This is the highest-level entry point for setting up a CMAP analysis
    environment from scratch.

    Parameters
    ----------
    source:
        Data source: ``"clue"`` for clue.io API, ``"legacy"`` for Broad
        Institute FTP, or ``"both"`` for both.
    output_dir:
        Directory for all downloaded and generated files.
    api_key:
        clue.io API key (required when *source* is ``"clue"`` or
        ``"both"``).
    n_permutations:
        Number of permutations for p-value null distributions.
    n_workers:
        Number of parallel workers for specificity computation.
    seed:
        Random seed for reproducibility.

    Returns
    -------
    Path
        The *output_dir* path with all installed data.

    Raises
    ------
    ValueError
        When *source* includes ``"clue"`` but no *api_key* is provided.
    """
    from pycmap.install import install  # noqa: PLC0415

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if source in ("clue", "both"):
        if api_key is None:
            raise ValueError(
                "An API key is required for clue.io. Register at https://clue.io to obtain one."
            )
        fetch_from_clue(api_key, output_dir)

    if source in ("legacy", "both"):
        fetch_legacy_data(output_dir)
        fetch_msigdb_gene_sets(output_dir)

    logger.info("Running CMAP installation pipeline...")
    install(
        data_dir=output_dir,
        n_permutations=n_permutations,
        n_workers=n_workers,
        seed=seed,
    )

    logger.info("Fetch and install complete. Data directory: %s", output_dir)
    return output_dir
