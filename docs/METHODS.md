# Methodology Documentation

This document describes the algorithms and methods used in the PLM vs O/E curve analysis.

## Overview

The PLM vs O/E curve analysis orders missense variants by their predicted pathogenicity (PLM score) and computes a sliding window observed/expected ratio along this ordering. This produces a curve showing how population constraint varies across the prediction score range.

## Key Concepts

### PLM Score (Protein Language Model Score)

PLM scores represent the predicted pathogenicity of each missense variant. Three models are available:

1. **Constraint Model**: Incorporates population constraint features
2. **Core Model**: Uses sequence-based features only (no structural data)
3. **Complete Model**: Full feature set including AlphaFold structural predictions

Higher scores indicate higher predicted pathogenicity.

### O/E Ratio (Observed/Expected)

The O/E ratio measures population constraint:
- **Observed**: Number of variants actually seen in the population
- **Expected**: Number of variants expected based on mutation rates

```
O/E = Σ(observed variants) / Σ(expected mutations)
```

- O/E < 1: Depleted (under negative selection)
- O/E ≈ 1: No constraint
- O/E > 1: Enriched (rare, could indicate positive selection or sequencing artifact)

---

## Algorithm: Sliding Window O/E by PLM Score

### Step 1: Data Preparation

1. Load variant-level data (RGC or gnomAD)
2. Filter to missense variants with valid PLM scores
3. Extract transcript identifiers from region field

### Step 2: Order Variants by PLM Score

For each transcript, variants are ordered by PLM score (ascending):
- Low PLM = benign predictions
- High PLM = pathogenic predictions

```
Position 1: Lowest PLM score (most benign)
Position 2: Next lowest
...
Position N: Highest PLM score (most pathogenic)
```

### Step 3: Build Neighbor Arrays

Using Hail scan operations, collect neighboring variants:

**Forward Pass** (descending index):
- Collect "following" neighbors (higher PLM scores)
- Store up to MAX_WINDOW=150 neighbors

**Backward Pass** (ascending index):
- Collect "previous" neighbors (lower PLM scores)
- Store up to MAX_WINDOW=150 neighbors

### Step 4: Compute Local O/E

For each variant position:

1. Combine arrays: `[previous] + [current] + [following]`
2. Sum observed counts: `window_obs = Σ(binary_presence_exomes_XX_XY)`
3. Sum expected rates: `window_exp = Σ(roulette_AR_MR_scaled_exomes_XX_XY)`
4. Compute ratio: `local_oe = window_obs / window_exp`

```python
# Pseudocode
all_values = previous_values + [current] + following_values
window_obs = sum(v.obs for v in all_values)
window_exp = sum(v.exp for v in all_values)
local_oe = window_obs / window_exp if window_exp >= MIN_EXPECTED else null
```

### Step 5: Edge Handling

At transcript boundaries:
- Fewer neighbors available on one side
- Window automatically expands to include available neighbors
- MIN_EXPECTED threshold ensures reliable O/E estimates

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MAX_WINDOW` | 150 | Maximum neighbors on each side (300 total) |
| `MIN_EXPECTED` | 10.0 | Minimum expected count for valid O/E |

---

## Curve Characterization

### Curve Categories

Genes are classified based on curve shape:

1. **Monotonic Decreasing** (plm_oe_corr < -0.5)
   - Higher PLM → Lower O/E
   - Expected pattern: pathogenic predictions are depleted

2. **Flat** (oe_range < 0.2)
   - O/E approximately constant across PLM range
   - Gene may not be under strong selection

3. **Monotonic Increasing** (plm_oe_corr > 0.3)
   - Higher PLM → Higher O/E
   - Rare, may indicate model issues

4. **Other**
   - Complex or variable patterns

### Key Metrics

| Metric | Description |
|--------|-------------|
| `plm_oe_corr` | Pearson correlation between PLM score and local O/E |
| `oe_range` | max(local_oe) - min(local_oe) |
| `gene_oe` | Gene-level O/E ratio |
| `n_variants` | Number of missense positions |

---

## ClinVar Validation

ClinVar pathogenic/likely pathogenic variants are overlaid on curves to validate:
- P/LP variants should cluster at high PLM scores
- P/LP variants should have lower local O/E

**Filtering Applied to ClinVar**:
1. Pathogenic or Likely_pathogenic clinical significance
2. Missense molecular consequence
3. Germline origin (exclude somatic-only)
4. No conflicting interpretations

---

## Mutation Rate Scaling

Expected mutation counts are derived from the Roulette mutation rate model:

1. **13-mer context**: Primary mutation rate lookup
2. **Hierarchical fallback**: 11-mer → 9-mer → 7-mer for rare contexts
3. **Scaling factor**: Adjusted to match observed synonymous variants

```
scaled_mutation_rate = roulette_13mer_rate × scaling_factor
```

The scaling factor is computed per-dataset (RGC vs gnomAD) using Poisson optimization to match observed synonymous variant counts.

---

## Data Source Comparison

### RGC vs gnomAD

| Aspect | RGC | gnomAD |
|--------|-----|--------|
| Sample size | ~400K exomes | 730K exomes + 76K genomes |
| Sex stratification | Combined XX_XY only | Full (XX, XY, XX_XY) |
| Population | Research cohort | Population reference |

For consistency, the PLM vs O/E analysis uses only `XX_XY` combined exome data from both sources.

---

## Output Files

### Per-Gene Statistics (`gene_stats.tsv`)

| Column | Description |
|--------|-------------|
| `transcript` | Transcript identifier |
| `n_variants` | Number of missense positions |
| `total_obs` | Total observed variants |
| `total_exp` | Total expected mutations |
| `gene_oe` | Gene-level O/E |
| `plm_oe_corr` | PLM-O/E correlation |
| `oe_range` | O/E range |
| `curve_category` | Classification |
| `n_clinvar` | ClinVar P/LP count |

### Curve Data (`curve_data.parquet`)

| Column | Description |
|--------|-------------|
| `transcript` | Transcript identifier |
| `plm_score` | PLM prediction score |
| `local_oe` | Sliding window O/E |
| `is_clinvar_pathogenic` | ClinVar P/LP status |
| `window_size` | Number of variants in window |

---

## References

- Roulette mutation rate model
- gnomAD constraint methodology
- ClinVar clinical interpretation guidelines
