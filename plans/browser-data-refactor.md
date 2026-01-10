# browser_data.ipynb Refactoring Plan

This document tracks refactoring changes for `browser_data.ipynb`.

**IMPORTANT:** Do NOT mark tasks complete until schema is verified using `metadata.json.gz` or parquet inspection.

---

## Quick Reference: Verification Commands

```bash
# Inspect Hail table schema (without loading into Hail)
zcat /path/to/table.ht/metadata.json.gz | python3 -m json.tool | grep -E "field_name_pattern"

# Inspect parquet schema
python3 -c "import polars as pl; df = pl.read_parquet('file.parquet', n_rows=0); print(df.columns)"

# Count columns matching pattern
python3 -c "import polars as pl; df = pl.read_parquet('file.parquet', n_rows=0); print(len([c for c in df.columns if 'pattern' in c]))"
```

---

## Current Schema Status

Based on merged HT schema inspection (2026-01-10):

| Category | Status | Notes |
|----------|--------|-------|
| RGC Summary/VIR/O/E | ✅ Present | No changes needed |
| Position fields | ✅ Present | No changes needed |
| Training struct | ✅ Present | No changes needed |
| dbNSFP struct | ✅ Present | No changes needed |
| Domains array | ✅ Present | No changes needed |
| Predictions (preds) | ⚠️ Uses tuples | Needs struct refactor |
| gnomAD Coverage | ⚠️ 2 of 12 | Need 10 more columns |
| gnomAD Constraint | ❌ Missing | ~471 columns needed |
| ClinVar | ⚠️ Wrong format | Replace struct with array |
| Variant Frequency | ❌ Missing | 3 columns needed |
| phyloP Scores | ❌ Missing | 2 columns needed |
| Stacked dbNSFP | ❌ Missing | 2 columns needed |

---

## Task 1: Convert Tuple Aggregations to Structs

**Status:** ❌ NOT COMPLETE
**Priority:** High

The schema shows `preds.Constraint` etc. use `array<tuple(str, float32, int32)>` not structs.

### 1.1 Constraint Predictions

**Current (tuples):**
```python
preds = hl.struct(
    Constraint = hl.agg.collect(hl.tuple([preds_ht.alleles[1], preds_ht.Constraint_1000_General_pred, preds_ht.Constraint_1000_General_n_pred])),
    ...
)
```

**Change to (structs):**
```python
preds = hl.struct(
    Constraint = hl.agg.collect(hl.struct(
        alt=preds_ht.alleles[1],
        pred=preds_ht.Constraint_1000_General_pred,
        n_pred=preds_ht.Constraint_1000_General_n_pred
    )),
    Core = hl.agg.collect(hl.struct(
        alt=preds_ht.alleles[1],
        pred=preds_ht.Core_1000_General_pred,
        n_pred=preds_ht.Core_1000_General_n_pred
    )),
    Complete = hl.agg.collect(hl.struct(
        alt=preds_ht.alleles[1],
        pred=preds_ht.Complete_1000_General_pred,
        n_pred=preds_ht.Complete_1000_General_n_pred
    ))
)
```

### Verification
```bash
# After export, check schema shows struct not tuple:
zcat output.ht/metadata.json.gz | python3 -m json.tool | grep -A5 "Constraint"
# Should show: "element_type": {"struct": {"fields": [{"name": "alt"...
# NOT: "element_type": {"tuple": ...
```

- [ ] Schema shows `array<struct{alt, pred, n_pred}>` not `array<tuple>`

---

## Task 2: Expand gnomAD Coverage Data

**Status:** ❌ NOT COMPLETE (have 2 of 12 columns)
**Priority:** High

### Data Sources

| Table | Path |
|-------|------|
| Exomes | `/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht` |
| Genomes | `/storage/zoghbi/data/sharing/hail_tables/gnomadV3_coverage_struct.ht` |

### Pre-Implementation: Inspect Source Schema
```bash
zcat /storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht/metadata.json.gz | python3 -m json.tool | grep -E "over_"
```

### Implementation

```python
coverage_thresholds = [10, 15, 20, 25, 30, 50]

gnomadV4_exomes_cov = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht')
gnomadV3_genomes_cov = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV3_coverage_struct.ht')

# Single join per table, then extract fields
base_ht = base_ht.annotate(
    _exomes_cov = gnomadV4_exomes_cov[base_ht.locus].gnomADV4_coverage,
    _genomes_cov = gnomadV3_genomes_cov[base_ht.locus].gnomADV3_coverage,
)

base_ht = base_ht.annotate(**{
    f'gnomad_exomes_over_{t}': base_ht._exomes_cov[f'over_{t}']
    for t in coverage_thresholds
}, **{
    f'gnomad_genomes_over_{t}': hl.float64(base_ht._genomes_cov[f'over_{t}'])
    for t in coverage_thresholds
})

base_ht = base_ht.drop('_exomes_cov', '_genomes_cov')
```

### Expected Columns (12)
- `gnomad_exomes_over_{10,15,20,25,30,50}`
- `gnomad_genomes_over_{10,15,20,25,30,50}`

### Verification
```python
# After export:
cols = [c for c in merged_ht.row if 'gnomad' in c and 'over_' in c]
print(f"Coverage columns: {len(cols)}")  # Should be 12
print(sorted(cols))
```

- [ ] Source table schema inspected
- [ ] 12 coverage columns present in output

---

## Task 3: Add ClinVar Variants Array

**Status:** ❌ NOT COMPLETE
**Priority:** High

Replace the old `clinvar` struct with a new `clinvar_variants` array.

### Current (wrong format)
```python
'clinvar': struct {
    clinvar_count: int64,
    clinvar_status_list: array<str>,
    clinvar_label_list: array<str>,
    clinvar_var_type_list: array<str>
}
```

### Target Format
```python
'clinvar_variants': array<struct {
    alt: str,            # Alternate allele
    significance: str,   # Pathogenic, Likely_pathogenic, VUS, etc.
    status: str,         # Review status (star rating)
    mol_csq: str,        # Molecular consequence
    variation_id: str    # ClinVar ID for linking
}>
```

### Implementation

```python
# Load ClinVar data (adjust path as needed)
clinvar_ht = hl.read_table('/path/to/clinvar.ht')

# Aggregate by locus and transcript
clinvar_agg = clinvar_ht.group_by('locus', 'transcript_id').aggregate(
    clinvar_variants = hl.agg.collect(hl.struct(
        alt = clinvar_ht.alleles[1],
        significance = clinvar_ht.clinical_significance,
        status = clinvar_ht.review_status,
        mol_csq = clinvar_ht.molecular_consequence,
        variation_id = hl.str(clinvar_ht.variation_id)
    ))
)

# Merge and drop old struct
merged_ht = merged_ht.annotate(
    clinvar_variants = clinvar_agg[merged_ht.locus, merged_ht.transcript_id].clinvar_variants
)
merged_ht = merged_ht.drop('clinvar')  # Remove deprecated struct
```

### Verification
```python
# After export:
assert 'clinvar_variants' in merged_ht.row
assert 'clinvar' not in merged_ht.row  # Old struct removed
merged_ht.select('clinvar_variants').show(3)
```

- [ ] Source ClinVar table schema inspected
- [ ] `clinvar_variants` column present
- [ ] Old `clinvar` struct removed

---

## Task 4: Add Variant Frequency Arrays

**Status:** ❌ NOT COMPLETE
**Priority:** High

### Target Schema
```python
'rgc_variants': array<struct{alt, af, ac, an, filters}>
'gnomad_exomes_variants': array<struct{alt, af, ac, an, filters}>
'gnomad_genomes_variants': array<struct{alt, af, ac, an, filters}>
```

### Data Sources

Identify where AF/AC/AN data comes from:
```bash
# Check RGC scaled table
zcat /storage/zoghbi/.../rgc_scaled.ht/metadata.json.gz | python3 -m json.tool | grep -E "(af|ac|an|AF|AC|AN)"

# Check gnomAD table
zcat /storage/zoghbi/.../gnomadV4_scaled.ht/metadata.json.gz | python3 -m json.tool | grep -E "(af|ac|an|AF|AC|AN)"
```

### Implementation

```python
# RGC variants (adjust field names after schema inspection)
rgc_ht = hl.read_table('/path/to/rgc_scaled.ht')
rgc_agg = rgc_ht.group_by('locus', 'transcript_id').aggregate(
    rgc_variants = hl.agg.collect(hl.struct(
        alt = rgc_ht.alleles[1],
        af = rgc_ht.info.AF,
        ac = rgc_ht.info.AC,
        an = rgc_ht.info.AN,
        filters = hl.if_else(hl.len(rgc_ht.filters) == 0, "PASS", hl.str(rgc_ht.filters))
    ))
)

# gnomAD exomes variants
gnomad_exomes_ht = hl.read_table('/path/to/gnomad_exomes.ht')
gnomad_exomes_agg = gnomad_exomes_ht.group_by('locus', 'transcript_id').aggregate(
    gnomad_exomes_variants = hl.agg.collect(hl.struct(
        alt = gnomad_exomes_ht.alleles[1],
        af = gnomad_exomes_ht.freq[0].AF,
        ac = gnomad_exomes_ht.freq[0].AC,
        an = gnomad_exomes_ht.freq[0].AN,
        filters = hl.if_else(hl.len(gnomad_exomes_ht.filters) == 0, "PASS", hl.str(gnomad_exomes_ht.filters))
    ))
)

# gnomAD genomes variants (similar pattern)
```

### Verification
```python
variant_cols = [c for c in merged_ht.row if c.endswith('_variants')]
print(f"Variant columns: {variant_cols}")  # Should be 3
merged_ht.select('rgc_variants').show(2)
```

- [ ] RGC source table schema inspected
- [ ] gnomAD source table schemas inspected
- [ ] `rgc_variants` column present
- [ ] `gnomad_exomes_variants` column present
- [ ] `gnomad_genomes_variants` column present

---

## Task 5: Add gnomAD Constraint Metrics

**Status:** ❌ NOT COMPLETE
**Priority:** High

### Data Source

**Table:** `/storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht`

### Pre-Implementation: Inspect Schema
```bash
zcat /storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht/metadata.json.gz | python3 -m json.tool | head -200
```

### Implementation

```python
gnomad_constraint_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht')

# Check key structure
print(f"Key: {list(gnomad_constraint_ht.key)}")
print(f"Row fields: {list(gnomad_constraint_ht.row)[:20]}...")

# Single join, then extract all constraint fields
base_ht = base_ht.annotate(
    _gnomad_constraint = gnomad_constraint_ht[base_ht.locus]
)

# Extract O/E columns dynamically
gnomad_constraint_cols = [col for col in gnomad_constraint_ht.row if '_oe' in col or 'vir' in col]
print(f"Found {len(gnomad_constraint_cols)} gnomAD constraint columns")

base_ht = base_ht.annotate(**{
    f'gnomad_{col}': base_ht._gnomad_constraint[col]
    for col in gnomad_constraint_cols
})

base_ht = base_ht.drop('_gnomad_constraint')
```

### Verification
```python
gnomad_oe_cols = [c for c in merged_ht.row if c.startswith('gnomad_') and '_oe' in c]
print(f"gnomAD O/E columns: {len(gnomad_oe_cols)}")  # Should be > 0
```

- [ ] Source table schema inspected
- [ ] Field naming pattern identified
- [ ] gnomAD constraint columns present in output (>100 expected)

---

## Task 6: Add phyloP Conservation Scores

**Status:** ❌ NOT COMPLETE
**Priority:** Medium

### Data Source

**Table:** `/storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht`

### Pre-Implementation: Inspect Schema
```bash
zcat /storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht/metadata.json.gz | python3 -m json.tool
```

### Implementation

```python
phylop_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht')

# Check field names (adjust based on schema inspection)
print(f"phyloP fields: {list(phylop_ht.row)}")

base_ht = base_ht.annotate(
    phylop_scores_447way = phylop_ht[base_ht.locus].phylop_447way,  # Adjust field name
    phylop_scores_100way = phylop_ht[base_ht.locus].phylop_100way,  # Adjust field name
)
```

### Expected Columns
- `phylop_scores_447way`
- `phylop_scores_100way`

### Verification
```python
phylop_cols = [c for c in merged_ht.row if c.startswith('phylop')]
print(f"phyloP columns: {phylop_cols}")  # Should be 2
```

- [ ] Source table schema inspected
- [ ] Exact field names identified
- [ ] `phylop_scores_447way` present
- [ ] `phylop_scores_100way` present

---

## Task 7: Add Stacked dbNSFP Tracks

**Status:** ❌ NOT COMPLETE
**Priority:** Medium

### Target Schema
```python
'AlphaMissense_stacked': array<struct{alt, score, percentile}>
'ESM1b_stacked': array<struct{alt, score, percentile}>
```

### Implementation

```python
# In the dbNSFP stacked aggregation cell
dbnsfp_stacked = hl.struct(
    AlphaMissense_stacked = hl.agg.collect(hl.struct(
        alt = ht.alleles[1],
        score = ht.AlphaMissense_am_pathogenicity,
        percentile = ht.AlphaMissense_am_pathogenicity_exome_perc
    )),
    ESM1b_stacked = hl.agg.collect(hl.struct(
        alt = ht.alleles[1],
        score = ht.ESM1b_score,
        percentile = ht.ESM1b_score_exome_perc
    ))
)
```

### Verification
```python
stacked_cols = [c for c in merged_ht.row if c.endswith('_stacked')]
print(f"Stacked columns: {stacked_cols}")  # Should be 2
```

- [ ] `AlphaMissense_stacked` present with struct schema
- [ ] `ESM1b_stacked` present with struct schema

---

## Task 8: Remove Deprecated Columns

**Status:** ❌ NOT COMPLETE
**Priority:** Low

### Columns to Remove

| Column | Reason |
|--------|--------|
| `rgc_any_max_af` | Max AF tracks removed from UI |
| `rgc_syn_max_af` | Max AF tracks removed from UI |
| `rgc_mis_max_af` | Max AF tracks removed from UI |
| `rgc_stop_gained_max_af` | Max AF tracks removed from UI |
| `Constraint` (top-level) | Duplicate of `preds.Constraint` |
| `Core` (top-level) | Duplicate of `preds.Core` |
| `Complete` (top-level) | Duplicate of `preds.Complete` |
| `clinvar` (struct) | Replaced by `clinvar_variants` |

### Implementation

```python
deprecated = [
    'rgc_any_max_af', 'rgc_syn_max_af', 'rgc_mis_max_af', 'rgc_stop_gained_max_af',
    'Constraint', 'Core', 'Complete',
    'clinvar',
]

for col in deprecated:
    if col in merged_ht.row:
        merged_ht = merged_ht.drop(col)
        print(f"Dropped: {col}")
```

### Verification
```python
for col in deprecated:
    assert col not in merged_ht.row, f"{col} should be removed"
print("All deprecated columns removed")
```

- [ ] All max_af columns removed
- [ ] Duplicate prediction columns removed
- [ ] Old clinvar struct removed

---

## Task 9: Export with Schema Verification

**Status:** ❌ NOT COMPLETE
**Priority:** Critical

This task ensures all columns are preserved during Hail→Parquet export.

### Implementation

```python
# BEFORE EXPORT: Verify schema
print("=" * 60)
print("PRE-EXPORT SCHEMA VERIFICATION")
print("=" * 60)

all_cols = list(merged_ht.row)
print(f"Total columns: {len(all_cols)}")

# Check each category
checks = {
    'gnomAD constraint': len([c for c in all_cols if c.startswith('gnomad_') and ('_oe' in c or 'vir' in c)]),
    'gnomAD coverage': len([c for c in all_cols if c.startswith('gnomad_') and 'over_' in c]),
    'Variant frequency': len([c for c in all_cols if c.endswith('_variants')]),
    'phyloP': len([c for c in all_cols if c.startswith('phylop')]),
    'Stacked dbNSFP': len([c for c in all_cols if c.endswith('_stacked')]),
}

for name, count in checks.items():
    status = '✓' if count > 0 else '✗'
    print(f"  {status} {name}: {count}")

# Force schema resolution with explicit select
key_fields = list(merged_ht.key)
non_key_cols = [c for c in all_cols if c not in key_fields]
merged_ht = merged_ht.select(*non_key_cols)

# Export
merged_repartitioned = merged_ht.repartition(1)
spark_df = merged_repartitioned.to_spark()
spark_df.write.mode('overwrite').parquet(output_file)

print(f"\nExported to: {output_file}")
```

### Post-Export Verification

```python
import polars as pl

df = pl.read_parquet(output_file, n_rows=0)
print(f"Parquet columns: {len(df.columns)}")

# Verify all expected columns survived export
expected = {
    'gnomAD coverage (12)': len([c for c in df.columns if 'gnomad' in c and 'over_' in c]) == 12,
    'clinvar_variants': 'clinvar_variants' in df.columns,
    'rgc_variants': 'rgc_variants' in df.columns,
    'gnomad_exomes_variants': 'gnomad_exomes_variants' in df.columns,
    'gnomad_genomes_variants': 'gnomad_genomes_variants' in df.columns,
    'phyloP (2)': len([c for c in df.columns if 'phylop' in c]) == 2,
    'Stacked (2)': len([c for c in df.columns if '_stacked' in c]) == 2,
    'No max_af': len([c for c in df.columns if 'max_af' in c]) == 0,
    'No old clinvar': 'clinvar' not in df.columns or 'clinvar_variants' in df.columns,
}

print("\nPost-export verification:")
for check, passed in expected.items():
    status = '✓' if passed else '✗ FAILED'
    print(f"  {status} {check}")
```

- [ ] Pre-export schema verification passed
- [ ] Parquet file created
- [ ] Post-export verification passed

---

## Task 10: Calculate Percentiles in Polars (Post-Export)

**Status:** ❌ NOT COMPLETE
**Priority:** High

Calculate percentiles in Polars using the `average` method (faster and handles ties correctly).

### 10.1 Exome-Wide Percentiles

For each score, calculate percentile across all non-null values:

**Scores to percentile:**
- `AlphaMissense_am_pathogenicity` → `AlphaMissense_am_pathogenicity_exome_perc`
- `ESM1b_score` → `ESM1b_score_exome_perc`
- `RGC_MTR_MTR` → `RGC_MTR_MTR_exome_perc`
- `Constraint_pred` → `Constraint_pred_exome_perc`
- `Core_pred` → `Core_pred_exome_perc`
- `Complete_pred` → `Complete_pred_exome_perc`

### 10.2 Cross-Score Comparable Percentiles

For "Comparable %iles" track section - percentile across positions where ALL scores are defined.

**Raw score columns needed in Hail export (include nulls):**
- `AlphaMissense_am_pathogenicity`
- `ESM1b_score`
- `RGC_MTR_MTR`
- `Constraint_pred`
- `Core_pred`
- `Complete_pred`

**Output columns (calculated in Polars):**
- `AlphaMissense_comparable_perc`
- `ESM1b_comparable_perc`
- `RGC_MTR_comparable_perc`
- `Constraint_comparable_perc`
- `Core_comparable_perc`
- `Complete_comparable_perc`

### Implementation (Polars script)

```python
import polars as pl

# Load parquet exported from Hail
df = pl.read_parquet('rgc_browser_data_merged.parquet')

# ============================================================
# 10.1 EXOME-WIDE PERCENTILES
# ============================================================
# Calculate percentile for each score independently (non-null values only)

exome_perc_configs = [
    ('AlphaMissense_am_pathogenicity', 'AlphaMissense_am_pathogenicity_exome_perc'),
    ('ESM1b_score', 'ESM1b_score_exome_perc'),
    ('RGC_MTR_MTR', 'RGC_MTR_MTR_exome_perc'),
    ('Constraint_pred', 'Constraint_pred_exome_perc'),
    ('Core_pred', 'Core_pred_exome_perc'),
    ('Complete_pred', 'Complete_pred_exome_perc'),
]

for score_col, perc_col in exome_perc_configs:
    if score_col in df.columns:
        df = df.with_columns(
            pl.col(score_col)
            .rank(method='average')
            .over(pl.col(score_col).is_not_null())  # Rank only non-null
            .truediv(pl.col(score_col).is_not_null().sum())
            .mul(100)
            .alias(perc_col)
        )
        print(f"✓ Calculated {perc_col}")

# ============================================================
# 10.2 CROSS-SCORE COMPARABLE PERCENTILES
# ============================================================
# Percentile across positions where ALL 6 scores are defined

comparable_scores = [
    'AlphaMissense_am_pathogenicity',
    'ESM1b_score',
    'RGC_MTR_MTR',
    'Constraint_pred',
    'Core_pred',
    'Complete_pred',
]

comparable_outputs = [
    'AlphaMissense_comparable_perc',
    'ESM1b_comparable_perc',
    'RGC_MTR_comparable_perc',
    'Constraint_comparable_perc',
    'Core_comparable_perc',
    'Complete_comparable_perc',
]

# Create mask for positions with all scores defined
all_defined_mask = pl.lit(True)
for col in comparable_scores:
    if col in df.columns:
        all_defined_mask = all_defined_mask & pl.col(col).is_not_null()

# Count positions in comparable set
n_comparable = df.filter(all_defined_mask).height
print(f"\nComparable set: {n_comparable:,} positions with all 6 scores defined")

# Calculate comparable percentiles
for score_col, perc_col in zip(comparable_scores, comparable_outputs):
    if score_col in df.columns:
        # Rank within the comparable subset only
        df = df.with_columns(
            pl.when(all_defined_mask)
            .then(
                pl.col(score_col)
                .rank(method='average')
                .over(all_defined_mask)
                .truediv(n_comparable)
                .mul(100)
            )
            .otherwise(None)
            .alias(perc_col)
        )
        print(f"✓ Calculated {perc_col}")

# ============================================================
# SAVE UPDATED PARQUET
# ============================================================
output_path = 'rgc_browser_data_merged_with_percentiles.parquet'
df.write_parquet(output_path)
print(f"\nSaved to: {output_path}")

# Verify new columns
new_cols = [c for c in df.columns if '_perc' in c]
print(f"Percentile columns: {len(new_cols)}")
for col in sorted(new_cols):
    non_null = df[col].drop_nulls().len()
    print(f"  {col}: {non_null:,} values")
```

### Track Tree Updates

Add new "Comparable %iles" section in `track_tree.py`:

```python
def build_comparable_percentiles_tree() -> Dict[str, Any]:
    """Build the Comparable Percentiles track section."""
    return {
        "label": "Comparable %iles",
        "children": [
            {"label": "AlphaMissense", "fieldId": "AlphaMissense_comparable_perc"},
            {"label": "ESM1b", "fieldId": "ESM1b_comparable_perc"},
            {"label": "MTR", "fieldId": "RGC_MTR_comparable_perc"},
            {"label": "Constraint", "fieldId": "Constraint_comparable_perc"},
            {"label": "Core", "fieldId": "Core_comparable_perc"},
            {"label": "Complete", "fieldId": "Complete_comparable_perc"},
        ]
    }
```

### Columns Required in Hail Export

Ensure these raw score columns are exported (even with nulls):

| Column | Description |
|--------|-------------|
| `AlphaMissense_am_pathogenicity` | Raw AlphaMissense score |
| `ESM1b_score` | Raw ESM1b score |
| `RGC_MTR_MTR` | Raw MTR score |
| `Constraint_pred` | Raw Constraint prediction |
| `Core_pred` | Raw Core prediction |
| `Complete_pred` | Raw Complete prediction |

### Verification

```python
# Check percentile columns exist and have reasonable distributions
perc_cols = [c for c in df.columns if c.endswith('_perc')]
print(f"Percentile columns: {len(perc_cols)}")

for col in perc_cols:
    stats = df.select(
        pl.col(col).min().alias('min'),
        pl.col(col).max().alias('max'),
        pl.col(col).mean().alias('mean'),
        pl.col(col).is_not_null().sum().alias('count')
    ).row(0)
    print(f"  {col}: min={stats[0]:.1f}, max={stats[1]:.1f}, mean={stats[2]:.1f}, n={stats[3]:,}")
```

- [ ] Exome-wide percentiles calculated (6 columns)
- [ ] Comparable percentiles calculated (6 columns)
- [ ] Track tree updated with "Comparable %iles" section
- [ ] Percentile distributions verified (0-100 range, mean ~50)

---

## Summary Checklist

Before marking this plan complete, ALL boxes must be checked:

### Data Tasks (Hail)
- [ ] Task 1: Predictions use structs (not tuples)
- [ ] Task 2: 12 gnomAD coverage columns
- [ ] Task 3: `clinvar_variants` array present
- [ ] Task 4: 3 variant frequency arrays present
- [ ] Task 5: gnomAD constraint columns present
- [ ] Task 6: 2 phyloP columns present
- [ ] Task 7: 2 stacked dbNSFP columns present
- [ ] Task 8: Deprecated columns removed
- [ ] Task 9: Export verification passed

### Data Tasks (Polars Post-Processing)
- [ ] Task 10: Exome-wide percentiles (6 `*_exome_perc` columns)
- [ ] Task 10: Comparable percentiles (6 `*_comparable_perc` columns)

### Downstream Updates (after data is ready)
- [ ] `scripts/preprocess_browser_data.py` handles new columns
- [ ] `browser/backend/track_tree.py` updated with "Comparable %iles" section
- [ ] Browser displays new tracks correctly

---

## Change Log

| Date | Change | Status |
|------|--------|--------|
| 2026-01-10 | Plan reset - all statuses corrected | ❌ Not started |
| 2026-01-10 | Merged preprocessing-fixes.md content | - |
| 2026-01-10 | Added verification commands | - |
| 2026-01-10 | Added Task 10: Polars percentile calculation | - |
| 2026-01-10 | Added "Comparable %iles" track section | - |
