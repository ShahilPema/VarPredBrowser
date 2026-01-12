# Plan: Upgrade ClinVar Track to gnomAD-Style with New Data Structure

## Goal
Modify the frontend to use the new variant-level ClinVar data structure and add gnomAD-style filtering/display.

## Data Structure (IMPLEMENTED 2026-01-09)

```python
'clinvar_variants': array<struct {
    alt: str,            # Alternate allele (e.g., "G")
    significance: str,   # Clinical significance (e.g., "Pathogenic", "Benign")
    status: str,         # Review status
    mol_csq: str,        # Molecular consequence for filtering
    variation_id: str    # ClinVar variation ID for linking
}>
```

**Note:** `clinvar_count` and old list-based columns removed.

## Current Frontend State
- File: `browser/frontend/index.html`
- Uses old `clinvar.clinvar_label_list` (array of strings)
- Has basic significance filtering (P/LP/VUS/LB/B)
- Missing: consequence filtering, star filtering, expand/collapse, variant details

---

## Implementation Plan

### Step 1: Update Data Field Reference
**Files:** `browser/backend/track_tree.py`, `browser/frontend/index.html`

- Update `track_tree.py`: fieldId from `clinvar.clinvar_label_list` → `clinvar_variants`
- Update `isClinVarStackedTrack()` function to check for `clinvar_variants`
- Update `drawClinVarStacked()` to read new structure

### Step 2: Add Significance Categories
Map `clinical_significance` to categories:
```javascript
const SIGNIFICANCE_GROUPS = {
    pathogenic: ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'],
    uncertain: ['Uncertain_significance', 'Conflicting_interpretations'],
    benign: ['Benign', 'Likely_benign', 'Benign/Likely_benign'],
    other: [] // everything else
};

function categorizeSignificance(sig) {
    if (!sig) return 'other';
    for (const [cat, terms] of Object.entries(SIGNIFICANCE_GROUPS)) {
        if (terms.some(t => sig.includes(t))) return cat;
    }
    return 'other';
}
```

### Step 3: Add Consequence Categories
Map `clinvar_mol_csq` to categories (gnomAD style):
```javascript
const CONSEQUENCE_CATEGORIES = {
    plof: ['frameshift_variant', 'stop_gained', 'splice_donor_variant',
           'splice_acceptor_variant', 'start_lost', 'stop_lost'],
    missense: ['missense_variant', 'inframe_insertion', 'inframe_deletion'],
    synonymous: ['synonymous_variant'],
    splice_region: ['splice_region_variant'],
    other: [] // everything else
};

function getConsequenceCategory(csq) {
    if (!csq) return 'other';
    for (const [cat, terms] of Object.entries(CONSEQUENCE_CATEGORIES)) {
        if (terms.includes(csq)) return cat;
    }
    return 'other';
}
```

### Step 4: Update Filter UI
ClinVar filter panel with two filter groups:

1. **Clinical Significance** - checkboxes with "only" buttons:
   - Pathogenic / Likely pathogenic
   - Uncertain significance
   - Benign / Likely benign
   - Other

2. **Consequence Type** - checkboxes with "only" buttons:
   - pLoF (frameshift, stop gained, splice)
   - Missense
   - Synonymous
   - Other

### Step 5: Update drawClinVarStacked() → drawClinVarVariants()
Rename and modify to:
1. Read `clinvar_variants` array
2. Filter variants by: significance category + consequence category
3. Draw individual variant markers (NO binned mode - expanded only)
4. Color variants by clinical significance (pathogenic=red, uncertain=orange, benign=blue, other=gray)
5. **Use different symbol shapes by consequence** (gnomAD style):
   - pLoF: Cross/X symbol (╳)
   - Missense: Triangle (▲)
   - Splice region: Diamond (◆)
   - Synonymous/Other: Circle (●)

### Step 6: Add Filtered Variant Count Display
Show count of visible variants after filtering:
- Display "Showing X of Y ClinVar variants" in the filter panel
- Update count when filters change

### Step 7: Update Tooltip
Show variant details on hover:
- Alternate allele (from `alt` field)
- Clinical significance (from `significance` field)
- Review status (from `status` field)
- Molecular consequence (from `mol_csq` field)
- **Clickable ClinVar link**: `https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/`

---

## Files to Modify

| File | Changes |
|------|---------|
| `browser/backend/track_tree.py` | Update fieldId to `clinvar_variants` |
| `browser/frontend/index.html` | Update rendering, filters, tooltip |

---

## Filter State Variables (to add)
```javascript
let clinvarSignificanceFilter = { pathogenic: true, uncertain: true, benign: true, other: true };
let clinvarConsequenceFilter = { plof: true, missense: true, synonymous: true, other: true };
```

---

## Summary of Changes

1. **Data Pipeline** (`browser_data.ipynb`): ✅ COMPLETE
   - Restructured ClinVar from list-based to per-variant structs
   - Fields: `alt`, `significance`, `status`, `mol_csq`, `variation_id`
   - Removed `clinvar_count` (not needed)

2. **Backend** (`track_tree.py`): PENDING
   - Change fieldId from `clinvar.clinvar_label_list` to `clinvar_variants`

3. **Frontend** (`index.html`): PENDING
   - Add helper functions: `getConsequenceCategory()`, `categorizeSignificance()`
   - Add filter state variables for significance and consequence
   - Update filter UI with significance + consequence checkboxes with "only" buttons
   - Rename `drawClinVarStacked()` → `drawClinVarVariants()` - draw individual variant markers (no binned mode)
   - **Use symbol shapes by consequence**: Cross (pLoF), Triangle (missense), Diamond (splice), Circle (other)
   - **Add variant count display**: "Showing X of Y ClinVar variants"
   - Update tooltip with alt, significance, status, consequence

---

## Features NOT Implemented (data not available)

These gnomAD features require data we don't have in our structure:
- HGVS notation (hgvsc, hgvsp)
- Submission details modal (submitter names, conditions)
- MedGen condition links
- gnomAD presence/frequency integration
- Star filtering (user requested to skip)

---

## Change Log

| Date | Change | Status |
|------|--------|--------|
| 2026-01-08 | Initial plan created | Pending |
| 2026-01-09 | Data restructuring implemented in browser_data.ipynb | ✅ Complete |
| 2026-01-09 | Field names: `alt`, `significance`, `status`, `mol_csq` | ✅ Complete |
| 2026-01-09 | `clinvar_count` removed (not needed) | ✅ Complete |
| - | Frontend implementation | Pending |
