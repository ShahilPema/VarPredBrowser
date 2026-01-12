# Add gnomAD Constraint Metrics

This plan describes integrating gnomAD constraint metrics as browser tracks, mirroring the existing RGC constraint implementation.

---

## Data Source

**Hail Table:** `/storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht`

---

## Pre-Implementation: Inspect Schema

Before implementing, inspect the hail table to understand available fields:

```python
import hail as hl

gnomad_constraint_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht')
gnomad_constraint_ht.describe()
gnomad_constraint_ht.show(5)
```

**Expected fields** (based on gnomAD constraint methodology):
- O/E ratios by consequence (missense, synonymous, LoF) and window size
- Expected variant counts
- Possibly LOEUF, pLI, or other gene-level constraint scores at locus level

---

## Implementation Overview

### 1. Data Pipeline (`browser_data.ipynb`)

Add gnomAD constraint annotation similar to RGC:

```python
# Load gnomAD constraint metrics
gnomad_constraint_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/constraint_metrics_by_locus_revft.ht')

# Annotate base table
# Field names TBD based on schema inspection
base_ht = base_ht.annotate(
    # O/E Ratios - example pattern (adjust field names after inspection)
    gnomad_mis_oe_3bp = gnomad_constraint_ht[base_ht.locus].mis_oe_3bp,
    gnomad_mis_oe_9bp = gnomad_constraint_ht[base_ht.locus].mis_oe_9bp,
    gnomad_mis_oe_21bp = gnomad_constraint_ht[base_ht.locus].mis_oe_21bp,
    gnomad_mis_oe_45bp = gnomad_constraint_ht[base_ht.locus].mis_oe_45bp,
    gnomad_mis_oe_93bp = gnomad_constraint_ht[base_ht.locus].mis_oe_93bp,

    gnomad_syn_oe_3bp = gnomad_constraint_ht[base_ht.locus].syn_oe_3bp,
    # ... etc for syn and any consequences

    # Expected values
    gnomad_mis_e_3bp = gnomad_constraint_ht[base_ht.locus].mis_e_3bp,
    # ... etc
)
```

### 2. Column Naming Convention - IMPLEMENTED

**Actual field naming pattern** (as of 2026-01-09):

```
gnomad_{consequence}_{cohort}_{sex}_{window}_oe
```

**Components:**
- **consequence**: `mis`, `syn`, `stop_gained`
- **cohort**: `exomes`, `genomes`, `combined`
- **sex**: `XX`, `XY`, or omitted for combined sex
- **window**: `3bp`, `9bp`, `21bp`, `45bp`, `93bp`

**Examples from schema:**
| Field Name | Type | Description |
|------------|------|-------------|
| `gnomad_mis_exomes_XX_21bp_oe` | float32 | Missense O/E, exomes, females, 21bp |
| `gnomad_syn_genomes_XY_45bp_oe` | float32 | Synonymous O/E, genomes, males, 45bp |
| `gnomad_stop_gained_exomes_3bp_oe` | float32 | Stop-gained O/E, exomes, combined sex, 3bp |
| `gnomad_mis_combined_XX_9bp_oe_weighted` | float64 | Missense weighted O/E, exomes+genomes, females, 9bp |

**Combined (weighted) fields:**
```
gnomad_{consequence}_combined_{sex}_{window}_oe_weighted
```

These represent sample-size-weighted averages of exome and genome O/E values.

**Comparison with RGC:**
| RGC Pattern | gnomAD Pattern |
|-------------|----------------|
| `rgc_mis_exomes_XX_XY_3bp_oe_af0epos00` | `gnomad_mis_exomes_3bp_oe` (combined sex) |
| `rgc_mis_exomes_XX_XY_21bp_oe_af1eneg06` | `gnomad_mis_exomes_21bp_oe` (no AF filter) |

**Key differences from RGC:**
- gnomAD has separate sex-stratified fields (`XX`, `XY`) in addition to combined
- gnomAD has both exomes and genomes (RGC has exomes only)
- gnomAD has weighted combined fields (`_combined_*_oe_weighted`)
- gnomAD has no AF filter suffix (constraint computed differently)

---

## 3. Track Tree Updates (`browser/backend/track_tree.py`)

### Option A: Add gnomAD as Sibling to RGC (Recommended)

Add a new top-level "gnomAD" section with O/E structure reflecting actual field names:

```python
def build_gnomad_oe_tree() -> Dict[str, Any]:
    """Build the gnomAD O/E Ratios tree section."""

    def window_node(consequence: str, cohort: str, sex: str, window: str) -> Dict[str, Any]:
        """Create node for a specific window."""
        # Build field name: gnomad_{consequence}_{cohort}_{sex}_{window}_oe
        if sex:
            field_id = f"gnomad_{consequence}_{cohort}_{sex}_{window}_oe"
        else:
            field_id = f"gnomad_{consequence}_{cohort}_{window}_oe"
        return {"label": window, "fieldId": field_id}

    def cohort_node(consequence: str, cohort: str, sex_options: list) -> Dict[str, Any]:
        """Create cohort node with sex/window hierarchy."""
        windows = ["3bp", "9bp", "21bp", "45bp", "93bp"]
        children = []
        for sex, sex_label in sex_options:
            sex_children = [window_node(consequence, cohort, sex, w) for w in windows]
            children.append({"label": sex_label, "children": sex_children})
        return {"label": cohort.capitalize(), "children": children}

    # Sex options: (suffix, label)
    sex_stratified = [("XX", "Female"), ("XY", "Male")]
    combined_sex = [("", "All")]

    consequences = [
        ("mis", "Missense"),
        ("syn", "Synonymous"),
        ("stop_gained", "Stop Gained"),
    ]

    return {
        "label": "O/E Ratios",
        "children": [
            {
                "label": csq_label,
                "children": [
                    cohort_node(csq, "exomes", sex_stratified + combined_sex),
                    cohort_node(csq, "genomes", sex_stratified + combined_sex),
                    # Combined weighted (exomes + genomes)
                    {
                        "label": "Combined (Weighted)",
                        "children": [
                            {
                                "label": sex_label,
                                "children": [
                                    {"label": w, "fieldId": f"gnomad_{csq}_combined_{sex}_{w}_oe_weighted"}
                                    for w in ["3bp", "9bp", "21bp", "45bp", "93bp"]
                                ],
                            }
                            for sex, sex_label in [("XX", "Female"), ("XY", "Male"), ("", "All")]
                        ],
                    },
                ],
            }
            for csq, csq_label in consequences
        ],
    }
```

### Option B: Add gnomAD Under Comparators

If gnomAD constraint is meant for comparison with RGC:

```python
# In build_track_tree(), add under "Comparators" > "Constraint":
{
    "label": "Constraint",
    "children": [
        {"label": "MTR (dbNSFP)", "fieldId": "dbnsfp.max_RGC_MTR_MTR"},
        {"label": "MTR Percentile", "fieldId": "dbnsfp.max_RGC_MTR_MTRpercentile_exome"},
        {"label": "CCR Residual Percentile", "fieldId": "dbnsfp.max_Non_Neuro_CCR_resid_pctile"},
        # New gnomAD entries
        {"label": "gnomAD O/E Mis 21bp", "fieldId": "gnomad_mis_21bp_oe"},
        {"label": "gnomAD O/E Syn 21bp", "fieldId": "gnomad_syn_21bp_oe"},
    ],
}
```

### Recommended: Option A (Full Parallel Structure)

Add gnomAD as a top-level section in `build_track_tree()`:

```python
def build_track_tree() -> Dict[str, Any]:
    track_tree = {
        "label": "Tracks",
        "children": [
            # ... existing children ...

            # Add gnomAD section after RGC
            {
                "label": "gnomAD",
                "children": [
                    build_gnomad_oe_tree(),
                    # Optionally add VIRs if available in gnomAD data
                ],
            },
        ],
    }
```

---

## 4. Preprocessing Updates (`scripts/preprocess_browser_data.py`)

The gnomAD columns will be automatically included if they follow the naming pattern and exist in the source parquet. Verify by checking:

```python
# In get_columns_to_keep()
gnomad_cols = [col for col in df.columns if col.startswith('gnomad_')]
print(f"Found {len(gnomad_cols)} gnomAD columns")
```

If columns need explicit inclusion, add to the column selection:

```python
# Add gnomAD constraint columns
gnomad_constraint_cols = [col for col in df.columns if col.startswith('gnomad_') and
                          any(x in col for x in ['_oe', '_e_', '_o_'])]
columns_to_keep.extend(gnomad_constraint_cols)
```

---

## 5. Frontend Updates

No frontend changes required if:
- Columns are numeric (Float64)
- Track tree entries have valid `fieldId` references
- Data flows through existing `/api/track-data` endpoint

The D3.js renderer will automatically handle new numeric tracks.

---

## Actual Columns (Implemented)

The gnomAD constraint columns are now in the schema. Here are the categories:

### O/E by Cohort and Sex (~90 columns)

**Pattern:** `gnomad_{consequence}_{cohort}_{sex}_{window}_oe` (float32)

| Consequence | Cohorts | Sex Options | Windows | Count |
|-------------|---------|-------------|---------|-------|
| mis | exomes, genomes | XX, XY | 3bp, 9bp, 21bp, 45bp, 93bp | 20 |
| syn | exomes, genomes | XX, XY | 3bp, 9bp, 21bp, 45bp, 93bp | 20 |
| stop_gained | exomes, genomes | XX, XY | 3bp, 9bp, 21bp, 45bp, 93bp | 20 |

**Also includes combined-sex fields** (no sex suffix):
- `gnomad_mis_exomes_3bp_oe`, `gnomad_syn_genomes_21bp_oe`, etc.

### Combined Weighted O/E (~45 columns)

**Pattern:** `gnomad_{consequence}_combined_{sex}_{window}_oe_weighted` (float64)

These combine exome and genome data with sample-size weighting.

| Consequence | Sex Options | Windows | Count |
|-------------|-------------|---------|-------|
| mis | XX, XY, (all) | 3bp, 9bp, 21bp, 45bp, 93bp | 15 |
| syn | XX, XY, (all) | 3bp, 9bp, 21bp, 45bp, 93bp | 15 |
| stop_gained | XX, XY, (all) | 3bp, 9bp, 21bp, 45bp, 93bp | 15 |

**Total: ~135 gnomAD O/E columns** (varies by exact sex stratification)

---

## Validation Checklist

- [ ] Inspect gnomAD constraint hail table schema
- [ ] Identify exact field names in source table
- [ ] Add annotation to `browser_data.ipynb`
- [ ] Verify columns in exported parquet
- [ ] Add `build_gnomad_oe_tree()` function to `track_tree.py`
- [ ] Add gnomAD section to track tree
- [ ] Test track display in browser
- [ ] Compare gnomAD vs RGC O/E values for sanity check

---

## Notes

- gnomAD constraint metrics are typically computed on gnomAD v4 exomes/genomes
- Field names in the source hail table may differ from expected pattern - adjust implementation accordingly
- Consider adding percentile columns if available (for percentile track display)
- If gnomAD has VIR-equivalent metrics, add `build_gnomad_vir_tree()` function

---

## Change Log

| Date | Change | Status |
|------|--------|--------|
| 2026-01-08 | Initial plan created | Pending |
| 2026-01-09 | Data fields implemented in browser_data.ipynb | ✅ Complete |
| 2026-01-09 | Updated plan with actual field naming convention | ✅ Complete |
| - | Track tree integration | Pending |
| - | Frontend testing | Pending |
