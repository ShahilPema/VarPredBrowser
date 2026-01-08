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

### 2. Column Naming Convention

Follow a parallel pattern to RGC columns:

| RGC Pattern | gnomAD Pattern |
|-------------|----------------|
| `rgc_mis_exomes_XX_XY_3bp_oe_af0epos00` | `gnomad_mis_3bp_oe` |
| `rgc_mis_exomes_XX_XY_3bp_e_af0epos00` | `gnomad_mis_3bp_e` |
| `rgc_syn_exomes_XX_XY_21bp_oe_af0epos00` | `gnomad_syn_21bp_oe` |

**Simplified naming** for gnomAD (no exomes/XX_XY/af suffix needed if single cohort):
- `gnomad_{consequence}_{window}_{metric}`
- consequence: `mis`, `syn`, `any`, `lof`
- window: `3bp`, `9bp`, `21bp`, `45bp`, `93bp`
- metric: `oe`, `e`, `o`

---

## 3. Track Tree Updates (`browser/backend/track_tree.py`)

### Option A: Add gnomAD as Sibling to RGC

Add a new top-level "gnomAD" section with parallel O/E structure:

```python
def build_gnomad_oe_tree() -> Dict[str, Any]:
    """Build the gnomAD O/E Ratios tree section."""

    def window_node(consequence: str, window: str) -> Dict[str, Any]:
        return {
            "label": window,
            "children": [
                {
                    "label": "O/E",
                    "fieldId": f"gnomad_{consequence}_{window}_oe",
                },
                {
                    "label": "Expected",
                    "fieldId": f"gnomad_{consequence}_{window}_e",
                },
            ],
        }

    windows = ["3bp", "9bp", "21bp", "45bp", "93bp"]

    return {
        "label": "O/E Ratios",
        "children": [
            {"label": "Missense", "children": [window_node("mis", w) for w in windows]},
            {"label": "Synonymous", "children": [window_node("syn", w) for w in windows]},
            {"label": "Any", "children": [window_node("any", w) for w in windows]},
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

## Expected New Columns

| Column | Type | Description |
|--------|------|-------------|
| `gnomad_mis_3bp_oe` | Float64 | Missense O/E ratio, 3bp window |
| `gnomad_mis_9bp_oe` | Float64 | Missense O/E ratio, 9bp window |
| `gnomad_mis_21bp_oe` | Float64 | Missense O/E ratio, 21bp window |
| `gnomad_mis_45bp_oe` | Float64 | Missense O/E ratio, 45bp window |
| `gnomad_mis_93bp_oe` | Float64 | Missense O/E ratio, 93bp window |
| `gnomad_syn_*_oe` | Float64 | Synonymous O/E ratios (5 windows) |
| `gnomad_any_*_oe` | Float64 | Any consequence O/E ratios (5 windows) |
| `gnomad_mis_*_e` | Float64 | Missense expected counts (5 windows) |
| `gnomad_syn_*_e` | Float64 | Synonymous expected counts (5 windows) |
| `gnomad_any_*_e` | Float64 | Any consequence expected counts (5 windows) |

**Total: ~30 new columns** (3 consequences x 5 windows x 2 metrics)

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
