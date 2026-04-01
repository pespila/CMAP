# pycmap-score

**Connectivity Map (CMAP) algorithm — KS-based drug-gene connectivity scoring**

A Python port of the Connectivity Map algorithm described in:

> Lamb J, Crawford ED, Peck D, et al. **The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease.** *Science* 313(5795):1929–1935, 2006. DOI: [10.1126/science.1132939](https://doi.org/10.1126/science.1132939)

This library computes connectivity scores between a user-supplied gene-expression signature (up/down-regulated gene lists) and a reference database of drug-treated expression profiles using a modified signed Kolmogorov-Smirnov statistic.

## Installation

```bash
pip install pycmap-score
```

With data-fetching support (for downloading data from clue.io):

```bash
pip install pycmap-score[fetch]
```

For development:

```bash
pip install pycmap-score[dev]
```

## Quick Start (with synthetic data)

No download required — generate a synthetic dataset and run the full pipeline immediately:

```python
from pycmap import compute_detailed_results, compute_permuted_results, synthetic_data, read_grp

# Generate synthetic data that includes gene IDs from the test files
all_genes = read_grp("test_files/readUp.grp") + read_grp("test_files/readDown.grp")
data = synthetic_data(extra_genes=all_genes)

# Run detailed results
det = compute_detailed_results("test_files/readUp.grp", "test_files/readDown.grp", data)
print(det.to_dataframe().head(10))

# Run permuted results
perm = compute_permuted_results(det, data, mode="by_name")
print(perm.to_dataframe())
```

You can also use the Alzheimer gene sets:

```python
all_genes = (
    read_grp("test_files/alzheimerUp.grp") +
    read_grp("test_files/alzheimerDown.grp")
)
data = synthetic_data(extra_genes=all_genes)
det = compute_detailed_results("test_files/alzheimerUp.grp", "test_files/alzheimerDown.grp", data)
print(det.to_dataframe().head(10))
```

### With real CMAP data

```python
from pycmap import cmap

# Requires real rank matrix + metadata (see Data Setup below)
detailed, by_name, by_cell = cmap("readUp.grp", "readDown.grp", data_dir="/path/to/data")
print(by_name.to_dataframe().loc["geldanamycin"])
```

## Synthetic Data API

`synthetic_data()` generates a realistic CMAP-shaped dataset for testing, demos, and development:

```python
from pycmap import synthetic_data

# Full-size (22,283 genes x 6,100 instances) — ~260 MB, takes ~2 seconds
data = synthetic_data()

# Smaller for quick experiments
data = synthetic_data(n_genes=5000, n_instances=200)

# Ensure your query gene IDs are present in the matrix
data = synthetic_data(extra_genes=["BRCA1", "TP53", "MYC", "CDKN1A"])

# Custom number of drugs (each gets multiple instances for permutation testing)
data = synthetic_data(n_drugs=8)

# Reproducible (same seed = same data)
data = synthetic_data(seed=123)
```

Parameters:
- `n_genes` — number of genes/rows (default: 22,283)
- `n_instances` — number of drug instances/columns (default: 6,100)
- `extra_genes` — gene IDs to guarantee are in the matrix (for `.grp` file compatibility)
- `n_drugs` — number of unique drug names (default: auto, up to 15)
- `seed` — random seed (default: 42)

## Lower-Level API

```python
from pycmap import (
    compute_detailed_results,
    compute_permuted_results,
    load_cmap_data,
    read_grp,
    ks_statistic,
)

# Load real data from disk
data = load_cmap_data("/path/to/data")

# Compute detailed results only
up_genes = read_grp("alzheimerUp.grp")
down_genes = read_grp("alzheimerDown.grp")
detailed = compute_detailed_results(up_genes, down_genes, data)

# Or use the vectorised KS statistic directly
import numpy as np
ranks = np.array([1, 5, 10, 20])
score = ks_statistic(ranks, n=100)
```

## Data Setup

The algorithm requires a pre-computed rank matrix and metadata files. You can use `synthetic_data()` for testing (see above) or real data:

### Option A: Download and install automatically

```python
from pycmap.fetch import fetch_and_install

# From legacy Broad Institute sources (no API key needed)
fetch_and_install(source="legacy", output_dir="data/")

# From clue.io (requires free registration at https://clue.io)
fetch_and_install(source="clue", output_dir="data/", api_key="YOUR_KEY")
```

### Option B: Use existing data

Place the following files in a directory:
- Rank matrix: `rankMatrix.npz` (or legacy `rankMatrix.txt.headn1` + `rankMatrix.txt.sed1d`)
- Metadata: `instID.data`, `instanceList.data`, `cellLineList.data`, `instanceNames.data`, `instanceCellNames.data`
- Pre-computed (optional): `pValue.npz`, `specificity.npz`

## Algorithm Overview

The CMAP algorithm scores how strongly a drug's expression profile connects to a query gene signature:

1. **Signed KS Statistic**: For each drug instance, compute a modified Kolmogorov-Smirnov statistic measuring whether query genes cluster at the top or bottom of the drug's ranked gene list.

2. **Score Combination**: Combine up- and down-regulation KS scores. If they have opposite signs (indicating a consistent connectivity pattern), the combined score is `k_up - k_down`. Same-sign scores yield 0.

3. **Normalisation**: Scores are normalised to [-1, 1] by dividing positive scores by the maximum and negative scores by the absolute minimum.

4. **Permuted Results**: Drugs are grouped by name (or name + cell line). For each group, an enrichment score, p-value (from 100,000 random permutations), and specificity (percentile rank against MSigDB gene sets) are computed.

For full mathematical details, see [docs/math.md](docs/math.md).

## Project Structure

```
src/pycmap/
├── __init__.py      # Public API: cmap() entry point
├── core.py          # Vectorised KS statistic and scoring
├── detailed.py      # Detailed results computation
├── permuted.py      # Permuted results (enrichment, p-value, specificity)
├── data.py          # File I/O (.grp, .gmt, rank matrix)
├── install.py       # Pre-computation pipeline
├── fetch.py         # Data fetching (clue.io, legacy sources)
└── _types.py        # CmapData, DetailedResults, PermutedResults
```

The original R implementation is preserved in [`legacy/`](legacy/).

## Development

```bash
# Install in dev mode
pip install -e ".[dev]"

# Run tests
pytest -v

# Lint
ruff check src/ tests/

# Format
ruff format src/ tests/

# Build
python -m build
```

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Citation

If you use this software, please cite the original CMAP paper:

```bibtex
@article{lamb2006connectivity,
  title={The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease},
  author={Lamb, Justin and Crawford, Emily D and Peck, David and Modell, Joshua W and Blat, Irene C and Wrobel, Matthew J and Lerner, Jim and Brunet, Jean-Philippe and Subramanian, Aravind and Ross, Kenneth N and others},
  journal={Science},
  volume={313},
  number={5795},
  pages={1929--1935},
  year={2006},
  publisher={American Association for the Advancement of Science},
  doi={10.1126/science.1132939}
}
```
