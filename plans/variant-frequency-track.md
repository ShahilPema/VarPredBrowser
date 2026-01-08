# Variant Frequency Track - Frontend Implementation Plan

## Overview

Implement gnomAD-style variant frequency visualization tracks in the GeneBrowser frontend, displaying allele frequency data from multiple cohorts with visual AF encoding and filter status toggle.

**Data dependency:** `browser-data-refactor.md` Section 3 (Allele Frequency Collection)

---

## Visual Encoding (gnomAD Style)

Based on `@gnomad/track-variants` from gnomad-browser-toolkit:

| Property | Encoding |
|----------|----------|
| X-position | Genomic position (via filtered index scale) |
| Ellipse height (ry) | Log-scaled AF: `scaleLog().domain([0.00001, 0.01]).range([3, 14])` |
| Ellipse width (rx) | Fixed, based on bar width |
| Fill color | Cohort: blue=exomes, green=genomes, purple=joint, orange=RGC |
| Fill style | Solid = PASS, Hollow/dashed = Filtered |
| Y-offset | Stacked for multiple variants at same position |

---

## Implementation

### File: `gosling_mvp/frontend/d3_viewer.html`

### 1. Add Constants

```javascript
// Allele frequency log scale (gnomAD style)
const AF_SCALE = d3.scaleLog()
    .domain([0.00001, 0.01])  // AF range
    .range([3, 14])            // Ellipse height (ry) pixels
    .clamp(true);

// Cohort colors
const COHORT_COLORS = {
    gnomad_exomes: '#4285F4',   // Blue
    gnomad_genomes: '#34A853',  // Green
    gnomad_joint: '#9C27B0',    // Purple
    rgc: '#FF6D00'              // Orange
};
```

### 2. Add Drawing Method to TrackCanvas Class

```javascript
drawVariantFrequency(positionsData, xScale, barWidth, cohort, showFiltered) {
    const ctx = this.ctx;
    const centerY = this.height / 2;
    const rx = Math.max(1.5, barWidth / 3);

    positionsData.forEach((d, i) => {
        const variants = d.allele_frequencies || [];
        const x = xScale(d.filtered_idx);

        variants.forEach((variant, vIdx) => {
            const cohortData = variant[cohort];
            if (!cohortData || cohortData.af === null || cohortData.af === undefined) return;

            // Check filter status
            const isFiltered = cohortData.filters && cohortData.filters.length > 0;
            if (isFiltered && !showFiltered) return;

            // Calculate ellipse height from AF (log scale)
            const ry = cohortData.af > 0 ? AF_SCALE(cohortData.af) : 3;

            // Offset multiple variants vertically
            const offsetY = (vIdx - (variants.length - 1) / 2) * (ry + 2);

            // Draw ellipse
            ctx.beginPath();
            ctx.ellipse(x + barWidth / 2, centerY + offsetY, rx, ry, 0, 0, 2 * Math.PI);

            if (isFiltered) {
                // Filtered: hollow with dashed stroke
                ctx.strokeStyle = COHORT_COLORS[cohort];
                ctx.setLineDash([2, 2]);
                ctx.lineWidth = 1.5;
                ctx.stroke();
                ctx.setLineDash([]);
            } else {
                // PASS: filled solid
                ctx.fillStyle = COHORT_COLORS[cohort];
                ctx.fill();
                ctx.strokeStyle = '#333';
                ctx.lineWidth = 0.5;
                ctx.stroke();
            }
        });
    });
}
```

### 3. Add Overlap Mode for Multiple Cohorts

```javascript
drawVariantFrequencyOverlap(positionsData, xScale, barWidth, cohorts, showFiltered) {
    // Draw rarest cohorts first (for visibility)
    const order = ['rgc', 'gnomad_exomes', 'gnomad_genomes', 'gnomad_joint'];
    const sortedCohorts = cohorts.sort((a, b) => order.indexOf(a) - order.indexOf(b));

    sortedCohorts.forEach(cohort => {
        this.ctx.globalAlpha = 0.75;
        this.drawVariantFrequency(positionsData, xScale, barWidth, cohort, showFiltered);
    });
    this.ctx.globalAlpha = 1.0;
}
```

### 4. Add Track Type Detection

```javascript
function isVariantFrequencyTrack(trackId) {
    return trackId.startsWith('af_') || trackId === 'allele_frequencies';
}

function getVariantFrequencyCohort(trackId) {
    if (trackId === 'af_all_overlap') return null;  // Overlap mode
    return trackId.replace('af_', '');
}
```

### 5. Update Rendering Logic

In the main render loop where tracks are drawn:

```javascript
if (isVariantFrequencyTrack(trackId)) {
    const cohort = getVariantFrequencyCohort(trackId);
    if (cohort) {
        canvas.drawVariantFrequency(positionsData, xScaleZoomed, barWidth, cohort, showFilteredVariants);
    } else {
        // Overlap mode
        canvas.drawVariantFrequencyOverlap(
            positionsData, xScaleZoomed, barWidth,
            ['gnomad_exomes', 'gnomad_genomes', 'gnomad_joint', 'rgc'],
            showFilteredVariants
        );
    }
}
```

### 6. Add Filter Toggle UI

```html
<div class="track-control">
    <label class="filter-toggle">
        <input type="checkbox" id="show-filtered-variants" checked onchange="toggleFilteredVariants(this.checked)">
        Show filtered variants
    </label>
</div>
```

```javascript
let showFilteredVariants = true;

function toggleFilteredVariants(show) {
    showFilteredVariants = show;
    renderChart();  // Re-render all tracks
}
```

### 7. Update Tooltip

```javascript
function formatVariantFrequencyTooltip(d, trackId) {
    const cohort = getVariantFrequencyCohort(trackId) || 'all';
    const variants = d.allele_frequencies || [];
    const cohorts = cohort === 'all'
        ? ['gnomad_exomes', 'gnomad_genomes', 'gnomad_joint', 'rgc']
        : [cohort];

    let html = `
        <div class="tooltip-title">Variant Frequencies</div>
        <div class="tooltip-row">Position: <span class="tooltip-value">${d.chrom}:${d.pos.toLocaleString()}</span></div>
        <div class="tooltip-row">Gene: <span class="tooltip-value">${d.gene_symbol || 'N/A'}</span></div>
    `;

    variants.forEach(v => {
        html += `<div style="margin-top: 6px; border-top: 1px solid #444; padding-top: 4px;">
            <strong>${v.alt}</strong>
        </div>`;

        cohorts.forEach(c => {
            const data = v[c];
            if (!data || data.af === null) return;

            const afDisplay = data.af < 0.0001
                ? data.af.toExponential(2)
                : data.af.toPrecision(3);
            const filterStatus = (data.filters && data.filters.length > 0)
                ? `<span style="color:#f44336">${data.filters.join(', ')}</span>`
                : '<span style="color:#4caf50">PASS</span>';
            const cohortLabel = c.replace('gnomad_', 'gnomAD ').replace('rgc', 'RGC');

            html += `
                <div class="tooltip-row" style="margin-left: 8px;">
                    <span style="color:${COHORT_COLORS[c]}">${cohortLabel}</span>
                </div>
                <div class="tooltip-row" style="margin-left: 16px;">
                    AF: <span class="tooltip-value">${afDisplay}</span> ${filterStatus}
                </div>
                <div class="tooltip-row" style="margin-left: 16px;">
                    AC/AN: <span class="tooltip-value">${data.ac.toLocaleString()} / ${data.an.toLocaleString()}</span>
                </div>
            `;
        });
    });

    return html;
}
```

---

## Backend Changes

### File: `gosling_mvp/backend.py`

Add to TRACK_TREE:

```python
{
    "label": "Variant Frequencies",
    "children": [
        {"label": "gnomAD Exomes", "fieldId": "af_gnomad_exomes", "type": "variant_frequency"},
        {"label": "gnomAD Genomes", "fieldId": "af_gnomad_genomes", "type": "variant_frequency"},
        {"label": "gnomAD Joint", "fieldId": "af_gnomad_joint", "type": "variant_frequency"},
        {"label": "RGC", "fieldId": "af_rgc", "type": "variant_frequency"},
        {"label": "All Cohorts (Overlapped)", "fieldId": "af_all_overlap", "type": "variant_frequency_overlap"},
    ]
}
```

---

## Files to Modify

| File | Changes |
|------|---------|
| `gosling_mvp/frontend/d3_viewer.html` | Add AF_SCALE, COHORT_COLORS, drawVariantFrequency(), tooltip, filter toggle |
| `gosling_mvp/backend.py` | Add TRACK_TREE entries for AF tracks |

---

## Verification

1. **Data available**: Confirm `allele_frequencies` field in parquet after running browser_data.ipynb
2. **Backend**: Check `/api/track-tree` returns new variant frequency tracks
3. **Frontend**:
   - Enable AF track for a gene
   - Verify ellipses render with height proportional to AF
   - Toggle filter checkbox - confirm filtered variants appear/disappear
   - Hover over variant - verify tooltip shows AF, AC/AN, filter status
   - Test overlap mode with all cohorts visible
4. **Edge cases**:
   - Position with no variants (should render nothing)
   - Position with multiple variants (should stack vertically)
   - Very rare variants (AF < 0.00001) should have minimum visible size

---

## Change Log

| Date | Change | Status |
|------|--------|--------|
| 2026-01-08 | Initial frontend implementation plan | Pending |
