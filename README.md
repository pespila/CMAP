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

## Quick Start

```python
from pycmap import cmap

# Run full analysis — returns detailed, by-name, and by-cell results
detailed, by_name, by_cell = cmap("readUp.grp", "readDown.grp", data_dir="/path/to/data")

# View top results as pandas DataFrames
print(detailed.to_dataframe().head(10))
print(by_name.to_dataframe().head(10))

# Query a specific drug
print(by_name.to_dataframe().loc["geldanamycin"])

# Query by drug + cell line
print(by_cell.to_dataframe().loc["tanespimycin - MCF7"])
```

## Lower-Level API

```python
from pycmap import (
    compute_detailed_results,
    compute_permuted_results,
    load_cmap_data,
    read_grp,
    ks_statistic,
)

# Load data
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

The algorithm requires a pre-computed rank matrix and metadata files. You can:

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
