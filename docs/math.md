# Mathematical Background: The Connectivity Map Algorithm

## Reference

Lamb J, Crawford ED, Peck D, et al. **The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease.** *Science* 313(5795):1929–1935, 2006. DOI: [10.1126/science.1132939](https://doi.org/10.1126/science.1132939)

---

## 1. Rank Matrix

The CMAP database consists of a **rank matrix** of dimensions *n* × *m*, where:
- *n* ≈ 22,000 genes (rows)
- *m* ≈ 6,100 drug treatment instances (columns)

Each column represents one drug treatment experiment. The values are **integer ranks** from 1 to *n*, representing the gene's expression level relative to all other genes under that treatment. Rank 1 = most up-regulated; rank *n* = most down-regulated.

## 2. Modified Signed Kolmogorov-Smirnov Statistic

The core of the CMAP algorithm is a **modified one-sided KS statistic** that differs from the standard two-sample KS test. It returns a **signed** value, encoding the direction of enrichment.

Given:
- **v** = sorted rank positions of query genes within a drug instance (1-indexed)
- **t** = number of query genes (length of v)
- **n** = total number of genes in the rank space

The statistic is:

```
j = 1, 2, ..., t

a = max_{j=1..t} (j/t − v_j/n)     [maximum positive deviation]
b = max_{j=1..t} (v_j/n − (j−1)/t)  [maximum negative deviation]

KS(v, n) = a    if a > b
         = −b   otherwise
```

**Interpretation:**
- **Positive KS**: Query genes are enriched toward the *top* of the ranked list (low rank values = high expression).
- **Negative KS**: Query genes are enriched toward the *bottom* of the ranked list.
- **Near zero**: No enrichment pattern.

This is computed separately for the up-regulated gene set (yielding *k_up*) and the down-regulated gene set (yielding *k_down*).

## 3. Combined Connectivity Score

For each drug instance, the up and down KS statistics are combined:

```
           ⎧ k_up − k_down    if sign(k_up) ≠ sign(k_down)
s_i =      ⎨
           ⎩ 0                 otherwise
```

A non-zero score requires **opposite signs**: up-regulated genes enriched in one direction, down-regulated genes enriched in the opposite. This indicates a consistent connectivity pattern between the query signature and the drug.

## 4. Score Normalisation

Raw scores are normalised to the interval [−1, 1]:

```
p = max(s_1, s_2, ..., s_m)     [maximum positive score]
q = min(s_1, s_2, ..., s_m)     [minimum negative score]

           ⎧  s_i / p      if s_i ≥ 0
ŝ_i =     ⎨
           ⎩ −(s_i / q)    if s_i < 0
```

The most positive score becomes +1.0, the most negative becomes −1.0.

## 5. Permuted Results

Detailed scores are aggregated by drug name (or drug name + cell line) to produce permuted results.

### 5.1 Enrichment Score

For each unique drug name, find the **row positions** of all its instances within the sorted detailed results. The enrichment score is the signed KS statistic applied to these positions:

```
Let V = sorted row positions (1-indexed) of instances matching drug name d
Let M = total number of instances

enrichment_d = KS(V, M)
```

This measures whether instances of drug *d* are clustered at the top (positive enrichment) or bottom (negative enrichment) of the ranked results.

### 5.2 P-value

The p-value is computed empirically from a null distribution of 100,000 random permutations:

```
For each unique frequency f (how many instances a drug has):
    Repeat 100,000 times:
        1. Draw f random positions from 1..M (without replacement)
        2. Sort them
        3. Compute |KS(positions, M)|
    Store the 100,000 absolute KS scores as the null distribution for frequency f

p_d = count(null_dist[f_d] ≥ |enrichment_d|) / 100,000
```

**Sentinel value**: p = 100 when:
- The drug appears only once (n = 1)
- The mean detailed score is exactly 0
- The non-null percentage is ≤ 50%

### 5.3 Non-Null Percentage

Measures the consistency of score direction across a drug's instances:

```
non_null_d = max(count(scores > 0), count(scores < 0)) / total × 100
```

Values close to 100% indicate consistent directionality.

### 5.4 Specificity

Specificity measures how specific the enrichment is to the query gene set compared to a reference panel of ~312 MSigDB gene sets.

For each drug, a **specificity matrix** is pre-computed: the enrichment score of every drug against every MSigDB gene set. Then:

```
For negative enrichment:
    specificity_d = count(enrichment_d ≥ negative_msig_scores) / count(negative_msig_scores)

For positive enrichment:
    specificity_d = 1 − count(enrichment_d ≥ positive_msig_scores) / count(positive_msig_scores)
```

**Lower specificity = more specific** (the drug's enrichment is more extreme compared to other pathways).

## 6. Output Matrices

### Detailed Results
| Column | Description |
|--------|-------------|
| InstID | Instance identifier |
| Score | Normalised connectivity score [−1, 1] |
| Up | KS statistic for up-regulated gene set |
| Down | KS statistic for down-regulated gene set |

Sorted by Score descending, then Up descending.

### Permuted Results
| Column | Description |
|--------|-------------|
| mean | Mean detailed score across instances |
| n | Number of instances |
| enrichment | KS enrichment score |
| p | Empirical p-value |
| specificity | Percentile rank vs MSigDB |
| % non-null | Score direction consistency |

Sorted by p ascending, then enrichment descending.

## 7. Why This Matters

The Connectivity Map enables **drug repurposing** and **mechanism of action discovery**:

- A disease gene signature queried against CMAP can identify drugs whose profiles **anti-correlate** with the disease (negative connectivity → potential therapeutic).
- Drugs with similar connectivity profiles may share **mechanisms of action**, even if structurally dissimilar.
- The approach is **unbiased** — it scans thousands of compounds simultaneously without requiring prior knowledge of drug targets.

The original 2006 paper demonstrated this by recovering known drug-disease connections (e.g., identifying HDAC inhibitors from a dexamethasone-resistance signature) and discovering novel connections subsequently validated experimentally.
