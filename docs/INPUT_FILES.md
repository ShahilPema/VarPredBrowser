# Input Files Documentation

This document describes all input data files used by the VarPredBrowser pipeline.

## Primary Data Sources

### 1. Variant-Level Scaled Data

These Hail tables contain variant-level data with mutation rates and observed variant counts.

#### rgc_scaled.ht

**Path**: `/storage/zoghbi/home/u235147/merged_vars/rgc_scaled.ht`
**Rows**: 817,452,381
**Source**: Regeneron Genetics Center (RGC)

**Key Columns**:
| Column | Type | Description |
|--------|------|-------------|
| `locus` | locus<GRCh38> | Genomic position |
| `alleles` | array<str> | Reference and alternate alleles |
| `region` | str | Transcript-position identifier (e.g., "ENST00000123456-100") |
| `most_deleterious_consequence_cds` | str | Consequence type (missense_variant, synonymous_variant, etc.) |
| `binary_presence_exomes_XX_XY` | int32 | Whether variant was observed (0 or 1) |
| `roulette_AR_MR_scaled_exomes_XX_XY` | float64 | Scaled mutation rate (expected) |
| `aa_pos` | int32 | Amino acid position |

**RGC-Specific Columns**:
| Column | Type | Description |
|--------|------|-------------|
| `rgc_call_rate` | float64 | Variant call rate in RGC data |
| `rgc_supmax` | float64 | Super maximum value from RGC |

---

#### gnomadV4_scaled.ht

**Path**: `/storage/zoghbi/home/u235147/merged_vars/gnomadV4_scaled.ht`
**Rows**: 817,452,381
**Source**: gnomAD v4 (Genome Aggregation Database)

**Key Columns** (same as rgc_scaled.ht):
| Column | Type | Description |
|--------|------|-------------|
| `locus` | locus<GRCh38> | Genomic position |
| `alleles` | array<str> | Reference and alternate alleles |
| `region` | str | Transcript-position identifier |
| `most_deleterious_consequence_cds` | str | Consequence type |
| `binary_presence_exomes_XX_XY` | int32 | Whether variant was observed |
| `roulette_AR_MR_scaled_exomes_XX_XY` | float64 | Scaled mutation rate |

**gnomAD-Specific Columns** (full sex stratification):
| Column | Type | Description |
|--------|------|-------------|
| `binary_presence_exomes_XX` | int32 | Observed in female exomes |
| `binary_presence_exomes_XY` | int32 | Observed in male exomes |
| `binary_presence_genomes_XX` | int32 | Observed in female genomes |
| `binary_presence_genomes_XY` | int32 | Observed in male genomes |
| `binary_presence_genomes_XX_XY` | int32 | Observed in combined genomes |
| `roulette_AR_MR_scaled_exomes_XX` | float64 | Scaled rate for female exomes |
| `roulette_AR_MR_scaled_exomes_XY` | float64 | Scaled rate for male exomes |
| `roulette_AR_MR_scaled_genomes_*` | float64 | Scaled rates for genomes |
| `gnomadV4_joint_freqs` | float64 | Joint allele frequency |

**Note**: For PLM vs O/E curves, we use only `XX_XY` combined exome columns for consistency.

---

### 2. Prediction Models

#### AOU_RGC_All_preds.ht

**Path**: `/local/Missense_Predictor_copy/Results/Inference/Predictions/AOU_RGC_All_preds.ht`
**Rows**: ~72 million (missense variants)
**Source**: ML models trained on pathogenicity data

**Prediction Columns**:
| Column | Description |
|--------|-------------|
| `Constraint_1000_General_aou_observed_neg4_pred` | **Constraint model** - uses population constraint features |
| `Core_1000_General_aou_observed_neg4_pred` | **Core model** - sequence-based features only |
| `Complete_1000_General_aou_observed_neg4_pred` | **Complete model** - full feature set with structural data |

**Score Interpretation**:
- Higher scores indicate higher predicted pathogenicity
- Scores are on a continuous scale
- Variants ordered by score for sliding window O/E calculation

---

### 3. ClinVar Data

#### cv_38_final.ht

**Path**: `/storage/zoghbi/data/sharing/hail_tables/cv_38_final.ht`
**Source**: ClinVar (NCBI)

**Key Fields** (nested under `clinvar_38` struct):
| Field | Description |
|-------|-------------|
| `CLNSIG` | Clinical significance (Pathogenic, Likely_pathogenic, etc.) |
| `MC` | Molecular consequence array |
| `ORIGIN` | Variant origin codes |

**Filtering Applied**:
1. Pathogenic or Likely pathogenic classification
2. Missense molecular consequence
3. Germline origin (excluding somatic-only)
4. No conflicting interpretations

---

### 4. Browser Data (Pre-merged)

#### rgc_browser_data_merged.parquet

**Path**: `/storage/zoghbi/home/u235147/merged_vars/rgc_browser_data_merged.parquet`
**Rows**: 33,387,608
**Columns**: 278

**Contents**:
- Merged locus-level data from all sources
- ClinVar annotations
- Training labels
- dbNSFP pathogenicity scores
- Domain annotations
- Constraint metrics (O/E, VIR)

**Used By**: Browser data pipeline and Gosling visualization

---

## Supporting Data Files

### Mutation Rate Tables (Roulette)

Hierarchical mutation rate lookup tables:
- `/storage/zoghbi/data/sharing/hail_tables/roulette_13mer.ht`
- `/storage/zoghbi/data/sharing/hail_tables/roulette_11mer.ht`
- `/storage/zoghbi/data/sharing/hail_tables/roulette_9mer.ht`
- `/storage/zoghbi/data/sharing/hail_tables/roulette_7mer.ht`

### Coverage Data

- `gnomadV4_exomes_coverage_struct.ht` - gnomAD v4 exome coverage
- `gnomadV3_coverage_struct.ht` - gnomAD v3 genome coverage

---

## Data Flow

```
rgc_scaled.ht / gnomadV4_scaled.ht
            ↓
    [Filter to missense]
            ↓
    AOU_RGC_All_preds.ht
            ↓
    [Join PLM scores]
            ↓
    cv_38_final.ht
            ↓
    [Annotate ClinVar status]
            ↓
    [Sliding window O/E calculation]
            ↓
    Per-gene PLM vs O/E curves
```

---

## File Sizes Summary

| File | Approximate Size |
|------|------------------|
| rgc_scaled.ht | ~100 GB |
| gnomadV4_scaled.ht | ~150 GB |
| AOU_RGC_All_preds.ht | ~50 GB |
| cv_38_final.ht | ~2 GB |
| browser_data_merged.parquet | ~21 GB |

**Total Storage Required**: ~300+ GB

---

## Access Requirements

These data files are stored on internal compute infrastructure and require:
- Access to `/storage/zoghbi/` mount
- Access to `/local/` directory for ML predictions
- Hail 0.2.134+ installed
- Sufficient memory for Spark processing (recommended: 40+ CPUs, 500+ GB RAM)
