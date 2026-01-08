# Plan: Add Allele Frequency Data to browser_data.ipynb

## Goal
Add gnomAD and Regeneron (RGC) allele frequency data (AC, AN, AF) to the browser data, structured as an array of per-allele frequency structs with filter status.

## Data Sources

| Source | Path | Key Fields |
|--------|------|------------|
| RGC | `/storage/zoghbi/data/sharing/hail_tables/no_anno/rgc.ht` | `info.ALL_AC`, `info.ALL_AF`, `info.ALL_AN`, `filters` |
| gnomAD Exomes | `/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes/gnomadV4exomes_snvs_sex.ht` | `gnomadV4_exomes.AC`, `.AN`, `.AF`, `.filters` |
| gnomAD Genomes | `/storage/zoghbi/data/sharing/hail_tables/gnomadV4_genomes/gnomadV4genomes_snvs.ht` | `gnomadV4_genomes.AC`, `.AN`, `.AF`, `.filters` |
| gnomAD Joint GrpMax | `/storage/zoghbi/data/sharing/hail_tables/no_anno/gnomad_joint_grpmax.ht` | `joint_grpmax`, `grpmax_anc` |

## Target Schema

```python
allele_frequencies: array<struct{
    alt: str,
    gnomad_exomes: struct{af: float64, ac: int64, an: int64, filtered: bool},
    gnomad_genomes: struct{af: float64, ac: int64, an: int64, filtered: bool},
    gnomad_joint_grpmax: float64,
    gnomad_grpmax_anc: str,
    rgc: struct{af: float64, ac: int64, an: int64, filtered: bool}
}>
```

## Implementation Plan

### Step 1: Load Frequency Tables

**File**: `notebooks/browser_data.ipynb`

Add a new cell after the existing data loading:

```python
# Load allele frequency tables
rgc_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/rgc.ht')
gnomad_exomes_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes/gnomadV4exomes_snvs_sex.ht')
gnomad_genomes_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_genomes/gnomadV4genomes_snvs.ht')
gnomad_joint_grpmax_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/gnomad_joint_grpmax.ht')
```

### Step 2: Join Frequency Data to Variant-Level Table

Join all frequency sources to the variant-level table (e.g., `all_sites` or equivalent):

```python
# Annotate with all frequency sources (keyed by locus, alleles)
freq_ht = all_sites.annotate(
    _rgc = rgc_ht[all_sites.key],
    _gnomad_exomes = gnomad_exomes_ht[all_sites.key].gnomadV4_exomes,
    _gnomad_genomes = gnomad_genomes_ht[all_sites.key].gnomadV4_genomes,
    _gnomad_joint = gnomad_joint_grpmax_ht[all_sites.key]
)
```

### Step 3: Build Frequency Struct in Select

Per user's preference, build the struct in `select()` before aggregation:

```python
# Select with frequency struct per allele
freq_ht = freq_ht.select(
    'locus',
    'alleles',
    'region',  # For grouping by transcript
    allele_freq = hl.struct(
        alt = freq_ht.alleles[1],
        gnomad_exomes = hl.struct(
            af = freq_ht._gnomad_exomes.AF,
            ac = freq_ht._gnomad_exomes.AC,
            an = freq_ht._gnomad_exomes.AN,
            filtered = hl.len(freq_ht._gnomad_exomes.filters) > 0
        ),
        gnomad_genomes = hl.struct(
            af = freq_ht._gnomad_genomes.AF,
            ac = freq_ht._gnomad_genomes.AC,
            an = freq_ht._gnomad_genomes.AN,
            filtered = hl.len(freq_ht._gnomad_genomes.filters) > 0
        ),
        gnomad_joint_grpmax = freq_ht._gnomad_joint.joint_grpmax,
        gnomad_grpmax_anc = freq_ht._gnomad_joint.grpmax_anc,
        rgc = hl.struct(
            af = freq_ht._rgc.info.ALL_AF,
            ac = freq_ht._rgc.info.ALL_AC,
            an = freq_ht._rgc.info.ALL_AN,
            filtered = hl.len(freq_ht._rgc.filters) > 0
        )
    )
)
```

### Step 4: Aggregate by Locus/Transcript

Group and collect the frequency structs:

```python
# Group by (locus, transcript) and collect frequency structs
freq_agg = freq_ht.group_by('locus', freq_ht.region.split('-')[0]).aggregate(
    allele_frequencies = hl.agg.collect(freq_ht.allele_freq)
)

freq_agg = freq_agg.checkpoint(f'{my_bucket}/tmp/allele_frequencies.ht', overwrite=True)
```

### Step 5: Merge into Browser Data

Join the aggregated frequency table to the merged browser data:

```python
# Key freq_agg appropriately
freq_agg = freq_agg.key_by('locus', 'transcript_id')

# Merge into the main browser table
merged = merged.annotate(
    allele_frequencies = freq_agg[merged.locus, merged.transcript_id].allele_frequencies
)
```

---

## Filter Handling Pattern

Based on existing code in `gnomadconstraint.ipynb`:

```python
# Filtered alleles are identified by non-empty filter set
# hl.len(all_sites.filters_exomes) > 0 indicates a filtered variant
# We keep these but mark them with filtered=True for display
```

This allows the browser to:
- Display all variants (filtered and passing)
- Visually distinguish filtered variants (e.g., grey out or add warning icon)
- Allow users to toggle filtered variant visibility

---

## Files to Modify

| Step | File | Changes |
|------|------|---------|
| 1-5 | `notebooks/browser_data.ipynb` | Add new cell for frequency data loading, struct building, and aggregation |
| - | `scripts/preprocess_browser_data.py` | Add `allele_frequencies` to exported columns (if needed) |

---

## Notes

- The `filtered` boolean is per-source (exomes, genomes, RGC may have different filter status)
- Missing data (variant not in source) will have null values for that source's struct
- This pattern matches existing stacked data (e.g., `AlphaMissense_stacked`, `ESM1b_stacked`)
