# browser_data.ipynb Refactoring Plan

This document tracks planned refactoring changes for `browser_data.ipynb`.

---

## 1. Convert Tuple Aggregations to Structs

**Status:** Pending
**Priority:** High
**Rationale:** Structs provide named field access, better Parquet/JSON compatibility, and self-documenting schemas.

### 1.1 dbNSFP Stacked Scores (cell id: cd7iyom5dcf)

**Current:**
```python
dbnsfp_stacked_ht = ht_stacked.group_by('locus', 'Ensembl_transcriptid').aggregate(
    dbnsfp_stacked = hl.struct(
        AlphaMissense = hl.agg.collect(hl.tuple([
            ht_stacked.alleles[1],
            ht_stacked.AlphaMissense_am_pathogenicity,
            ht_stacked.AlphaMissense_am_pathogenicity_exome_perc
        ])),
        ESM1b = hl.agg.collect(hl.tuple([
            ht_stacked.alleles[1],
            ht_stacked.ESM1b_score,
            ht_stacked.ESM1b_score_exome_perc
        ]))
    )
)
```

**Refactor to:**
```python
dbnsfp_stacked_ht = ht_stacked.group_by('locus', 'Ensembl_transcriptid').aggregate(
    dbnsfp_stacked = hl.struct(
        AlphaMissense = hl.agg.collect(hl.struct(
            alt=ht_stacked.alleles[1],
            score=ht_stacked.AlphaMissense_am_pathogenicity,
            percentile=ht_stacked.AlphaMissense_am_pathogenicity_exome_perc
        )),
        ESM1b = hl.agg.collect(hl.struct(
            alt=ht_stacked.alleles[1],
            score=ht_stacked.ESM1b_score,
            percentile=ht_stacked.ESM1b_score_exome_perc
        ))
    )
)
```

### 1.2 Constraint Predictions (cell id: knphdflxbh)

**Current:**
```python
preds_agg = preds_ht.group_by('Ensembl_transcriptid', 'locus').aggregate(
    preds = hl.struct(
        Constraint = hl.agg.collect(hl.tuple([preds_ht.alleles[1], preds_ht.Constraint_1000_General_pred, preds_ht.Constraint_1000_General_n_pred])),
        Core = hl.agg.collect(hl.tuple([preds_ht.alleles[1], preds_ht.Core_1000_General_pred, preds_ht.Core_1000_General_n_pred])),
        Complete = hl.agg.collect(hl.tuple([preds_ht.alleles[1], preds_ht.Complete_1000_General_pred, preds_ht.Complete_1000_General_n_pred]))
    )
)
```

**Refactor to:**
```python
preds_agg = preds_ht.group_by('Ensembl_transcriptid', 'locus').aggregate(
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
)
```

### 1.3 Variant Consequences (cell id: du9ru4eihwo)

**Current:**
```python
csq_agg = csq_ht.group_by('transcript_id', 'locus').aggregate(
    variant_consequences = hl.agg.collect(hl.tuple([
        csq_ht.alleles[1],
        csq_ht.csq_category
    ]))
)
```

**Refactor to:**
```python
csq_agg = csq_ht.group_by('transcript_id', 'locus').aggregate(
    variant_consequences = hl.agg.collect(hl.struct(
        alt=csq_ht.alleles[1],
        consequence=csq_ht.csq_category
    ))
)
```

---

## Naming Conventions for Collected Structs

Use consistent field names across all variant-level collected structs:

| Field | Description | Example |
|-------|-------------|---------|
| `alt` | Alternate allele (from `alleles[1]`) | `"T"` |
| `pred` | Prediction value (0-1 probability) | `0.85` |
| `n_pred` | Number of predictions / ensemble count | `3` |
| `score` | Raw score value | `-2.5` |
| `percentile` | Percentile/rank value (0-100) | `95.2` |
| `consequence` | Consequence category | `"missense"` |

---

## Schema Changes Summary

After refactoring, the schema will change from tuples to named structs:

| Field | Before | After |
|-------|--------|-------|
| `preds.Constraint` | `array<tuple(str, float32, int32)>` | `array<struct{alt: str, pred: float32, n_pred: int32}>` |
| `preds.Core` | `array<tuple(str, float32, int32)>` | `array<struct{alt: str, pred: float32, n_pred: int32}>` |
| `preds.Complete` | `array<tuple(str, float32, int32)>` | `array<struct{alt: str, pred: float32, n_pred: int32}>` |
| `dbnsfp_stacked.AlphaMissense` | `array<tuple(str, float64, float64)>` | `array<struct{alt: str, score: float64, percentile: float64}>` |
| `dbnsfp_stacked.ESM1b` | `array<tuple(str, float64, float64)>` | `array<struct{alt: str, score: float64, percentile: float64}>` |
| `variant_consequences` | `array<tuple(str, str)>` | `array<struct{alt: str, consequence: str}>` |

---

## Downstream Impact

Files that may need updates after refactoring:

1. **gosling_mvp/preprocess_mis_all.py** - Check for tuple indexing (`x[0]`, `x[1]`)
2. **Frontend JavaScript** - Update any code parsing the Parquet/JSON output
3. **API endpoints** - Verify response serialization

### Migration Pattern

```python
# Before (tuple access)
for variant in row.preds.Constraint:
    alt_allele = variant[0]
    prediction = variant[1]
    n_pred = variant[2]

# After (struct access)
for variant in row.preds.Constraint:
    alt_allele = variant.alt
    prediction = variant.pred
    n_pred = variant.n_pred
```

---

## Validation Checklist

After applying changes:

- [ ] Run `merged.describe()` to verify struct field names
- [ ] Export chr2 test parquet and validate with Polars
- [ ] Check gosling_mvp preprocessing script compatibility
- [ ] Verify frontend displays data correctly

---

## 2. Expand gnomAD Coverage Data

**Status:** Pending
**Priority:** Medium
**Rationale:** Support overlayed exome/genome coverage tracks with multiple metrics for visualization.

### Data Sources

Located at `/storage/zoghbi/data/sharing/hail_tables/`:

| Table | Key | Struct Name |
|-------|-----|-------------|
| `gnomadV4_exomes_coverage_struct.ht` | locus (GRCh38) | `gnomADV4_coverage` |
| `gnomadV3_coverage_struct.ht` | locus (GRCh38) | `gnomADV3_coverage` |

Available fields in coverage structs:
- `mean` (Float64)
- `median_approx` (Int32)
- `total_DP` (Int64)
- `over_1`, `over_5`, `over_10`, `over_15`, `over_20`, `over_25`, `over_30`, `over_50`, `over_100`

### Current Implementation (cell id: f76f3449)

```python
gnomadV4_exomes_coverage = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht')
gnomadV4_genomes_coverage = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV3_coverage_struct.ht')

base_ht = base_ht.annotate(
    gnomad_exomes_over_20 = gnomadV4_exomes_coverage[base_ht.locus].gnomADV4_coverage.over_20,
    gnomad_genomes_over_20 = gnomadV4_genomes_coverage[base_ht.locus].gnomADV3_coverage.over_20
)
```

### Refactor to:

```python
gnomadV4_exomes_coverage = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht')
gnomadV4_genomes_coverage = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV3_coverage_struct.ht')

base_ht = base_ht.annotate(
    # Exome coverage
    gnomad_exomes_mean = gnomadV4_exomes_coverage[base_ht.locus].gnomADV4_coverage.mean,
    gnomad_exomes_median = hl.float64(gnomadV4_exomes_coverage[base_ht.locus].gnomADV4_coverage.median_approx),
    gnomad_exomes_over_20 = gnomadV4_exomes_coverage[base_ht.locus].gnomADV4_coverage.over_20,
    gnomad_exomes_over_30 = gnomadV4_exomes_coverage[base_ht.locus].gnomADV4_coverage.over_30,
    # Genome coverage
    gnomad_genomes_mean = gnomadV4_genomes_coverage[base_ht.locus].gnomADV3_coverage.mean,
    gnomad_genomes_median = hl.float64(gnomadV4_genomes_coverage[base_ht.locus].gnomADV3_coverage.median_approx),
    gnomad_genomes_over_20 = hl.float64(gnomadV4_genomes_coverage[base_ht.locus].gnomADV3_coverage.over_20),
    gnomad_genomes_over_30 = hl.float64(gnomadV4_genomes_coverage[base_ht.locus].gnomADV3_coverage.over_30),
)
```

### New Columns

| Column | Type | Description |
|--------|------|-------------|
| `gnomad_exomes_mean` | Float64 | Mean exome coverage depth |
| `gnomad_exomes_median` | Float64 | Median exome coverage depth |
| `gnomad_exomes_over_20` | Float64 | Fraction with ≥20x exome coverage |
| `gnomad_exomes_over_30` | Float64 | Fraction with ≥30x exome coverage |
| `gnomad_genomes_mean` | Float64 | Mean genome coverage depth |
| `gnomad_genomes_median` | Float64 | Median genome coverage depth |
| `gnomad_genomes_over_20` | Float64 | Fraction with ≥20x genome coverage |
| `gnomad_genomes_over_30` | Float64 | Fraction with ≥30x genome coverage |

### Notes

- Cast `median_approx` (Int32) to Float64 for consistency
- Cast genome `over_X` fields from Float32 to Float64 for consistency
- Tables are already loaded in notebook - just need to extract more fields

### Downstream Impact

1. Re-run notebook to regenerate `rgc_browser_data_merged.parquet`
2. Update `scripts/preprocess_browser_data.py` to include new columns in axis tables
3. Backend/frontend changes per `plans/coverage-track.md`

---

## 3. Add Allele Frequency Collection

**Status:** Pending
**Priority:** High
**Rationale:** Support variant frequency tracks with AF/AC/AN display and filter status for multiple cohorts.

### Data Sources

Allele frequency data from four cohorts:
- **gnomAD v4 Exomes**
- **gnomAD v4 Genomes**
- **gnomAD v4 Joint** (exome + genome combined)
- **RGC** (Regeneron Genetics Center)

### Implementation

Add new aggregation cell to collect per-variant AF data:

```python
# ALLELE FREQUENCIES (Collect AF/AC/AN per variant for each cohort)
# Group by locus and transcript, collect structs with cohort-specific frequency data

af_agg = af_ht.group_by('locus', 'Ensembl_transcriptid').aggregate(
    allele_frequencies = hl.agg.collect(hl.struct(
        alt=af_ht.alleles[1],
        gnomad_exomes=hl.struct(
            af=af_ht.gnomad_exomes_af,
            ac=af_ht.gnomad_exomes_ac,
            an=af_ht.gnomad_exomes_an,
            filters=af_ht.gnomad_exomes_filters  # array<str>, empty = PASS
        ),
        gnomad_genomes=hl.struct(
            af=af_ht.gnomad_genomes_af,
            ac=af_ht.gnomad_genomes_ac,
            an=af_ht.gnomad_genomes_an,
            filters=af_ht.gnomad_genomes_filters
        ),
        gnomad_joint=hl.struct(
            af=af_ht.gnomad_joint_af,
            ac=af_ht.gnomad_joint_ac,
            an=af_ht.gnomad_joint_an,
            filters=af_ht.gnomad_joint_filters
        ),
        rgc=hl.struct(
            af=af_ht.rgc_af,
            ac=af_ht.rgc_ac,
            an=af_ht.rgc_an,
            filters=af_ht.rgc_filters
        )
    ))
)

af_agg = af_agg.checkpoint(f'{my_bucket}/tmp/allele_frequencies.ht', overwrite=True)
```

### Schema

```
allele_frequencies: array<struct{
    alt: str,
    gnomad_exomes: struct{af: float64, ac: int64, an: int64, filters: array<str>},
    gnomad_genomes: struct{af: float64, ac: int64, an: int64, filters: array<str>},
    gnomad_joint: struct{af: float64, ac: int64, an: int64, filters: array<str>},
    rgc: struct{af: float64, ac: int64, an: int64, filters: array<str>}
}>
```

### Field Definitions

| Field | Type | Description |
|-------|------|-------------|
| `alt` | str | Alternate allele |
| `af` | float64 | Allele frequency (AC/AN) |
| `ac` | int64 | Allele count |
| `an` | int64 | Allele number (total alleles) |
| `filters` | array<str> | Filter flags (empty array = PASS) |

### Filter Status Convention

- `filters = []` (empty array) → Variant passed all filters (PASS)
- `filters = ["AC0", "RF"]` → Variant failed listed filters

### Merge into Final Table

```python
# Merge allele frequencies (by locus, transcript)
af_ht = hl.read_table(f'{my_bucket}/tmp/allele_frequencies.ht')
af_ht = af_ht.key_by('locus', 'Ensembl_transcriptid')
merged = merged.annotate(
    allele_frequencies = af_ht[merged.locus, merged.transcript_id].allele_frequencies
)
```

### Notes

- Filter field uses array<str> to support multiple filter flags
- Empty filters array indicates PASS status (consistent with VCF convention)
- All cohorts use same struct format for consistency
- Missing data (variant not in cohort) will have null values

### Downstream Impact

1. Frontend: New variant frequency track type (see `plans/variant-frequency-track.md`)
2. Backend: Add track tree entries for AF tracks
3. Tooltip: Display AF in scientific notation, show AC/AN, filter status

---

## 4. Add phyloP Conservation Scores

**Status:** Pending
**Priority:** High
**Rationale:** phyloP tracks are defined in track_tree.py but no data is currently being loaded. Conservation scores are essential for variant interpretation.

### Data Source

Located at `/storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht`

### Current Issue

The browser defines three phyloP tracks in `track_tree.py` (lines 234-239):
- `phylop_scores_447way`
- `phylop_scores_100way`
- `phylop17way`

However, these columns are not present in the browser data parquet files, causing empty tracks.

### Implementation

Add phyloP annotation to base table:

```python
# Load phyloP scores
phylop_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht')

# Annotate base table with phyloP scores
# (Field names TBD - need to inspect table schema)
base_ht = base_ht.annotate(
    phylop_scores_447way = phylop_ht[base_ht.locus].phylop_447way,
    phylop_scores_100way = phylop_ht[base_ht.locus].phylop_100way,
    phylop17way = phylop_ht[base_ht.locus].phylop_17way
)
```

### Pre-Implementation Steps

1. **Inspect table schema** to get exact field names:
   ```python
   phylop_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/phyloPscores_hg38_final.ht')
   phylop_ht.describe()
   ```

2. **Verify key compatibility** - confirm table is keyed by locus (GRCh38)

3. **Check data types** - phyloP scores are typically Float32/Float64

### Expected Columns

| Column | Type | Description |
|--------|------|-------------|
| `phylop_scores_447way` | Float64 | 447-way vertebrate conservation |
| `phylop_scores_100way` | Float64 | 100-way vertebrate conservation |
| `phylop17way` | Float64 | 17-way primate conservation |

### Downstream Impact

1. **preprocess_browser_data.py** - Columns will automatically be included if present in source parquet
2. **Backend** - No changes needed (track definitions already exist)
3. **Frontend** - No changes needed (tracks already configured)

### Validation

- [ ] Inspect phyloP hail table schema
- [ ] Add annotation to browser_data.ipynb
- [ ] Verify columns present in exported parquet
- [ ] Confirm tracks display in browser

---

## Future Additions

_Add additional refactoring tasks below as needed._

---

## Change Log

| Date | Change | Status |
|------|--------|--------|
| 2026-01-07 | Initial plan: tuple → struct refactoring | Pending |
| 2026-01-07 | Add gnomAD coverage expansion | Pending |
| 2026-01-08 | Add allele frequency collection | Pending |
| 2026-01-08 | Add phyloP conservation scores | Pending |
