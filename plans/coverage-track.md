# Plan: Add gnomAD-style Coverage Tracks to VarPredBrowser

## Goal
Add coverage tracks with overlayed exome and genome distributions, similar to gnomAD browser's implementation.

## gnomAD Coverage Track Summary
- **Data**: Per-position metrics (mean, median, over_1, over_5, over_10, over_15, over_20, over_25, over_30, over_50, over_100)
- **Visualization**: Area charts with semi-transparent overlays (exome: steel blue @ 0.7 opacity, genome: green @ 0.5 opacity)
- **Adaptive rendering**: Bar charts for <100bp, area charts for larger regions

## Existing VarPredBrowser Data
From `config/paths.yaml`:
```yaml
gnomad_exomes_coverage: "/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes_coverage_struct.ht"
gnomad_genomes_coverage: "/storage/zoghbi/data/sharing/hail_tables/gnomadV3_coverage_struct.ht"
```

Current columns already in axis tables:
- `gnomad_exomes_over_20`
- `gnomad_genomes_over_20`

---

## Implementation Plan

### Phase 1: Data Preparation (Hail notebook)

**File**: `notebooks/browser_data.ipynb`

Add coverage metrics to the merged browser data:

```python
# Load gnomAD coverage tables
exome_cov = hl.read_table(paths['gnomad_exomes_coverage'])
genome_cov = hl.read_table(paths['gnomad_genomes_coverage'])

# Join coverage to main data by locus
browser_data = browser_data.annotate(
    coverage = hl.struct(
        exome = hl.struct(
            mean = exome_cov[browser_data.locus].mean,
            median = exome_cov[browser_data.locus].median,
            over_20 = exome_cov[browser_data.locus].over_20,
            over_30 = exome_cov[browser_data.locus].over_30
        ),
        genome = hl.struct(
            mean = genome_cov[browser_data.locus].mean,
            median = genome_cov[browser_data.locus].median,
            over_20 = genome_cov[browser_data.locus].over_20,
            over_30 = genome_cov[browser_data.locus].over_30
        )
    )
)
```

**Alternative**: If storing as separate columns is preferred:
```python
browser_data = browser_data.annotate(
    gnomad_exomes_mean = exome_cov[browser_data.locus].mean,
    gnomad_exomes_median = exome_cov[browser_data.locus].median,
    gnomad_genomes_mean = genome_cov[browser_data.locus].mean,
    gnomad_genomes_median = genome_cov[browser_data.locus].median,
    # over_20 columns already exist
)
```

### Phase 2: Preprocessing Script Update

**File**: `scripts/preprocess_browser_data.py`

Add coverage columns to the selection:

```python
COVERAGE_COLUMNS = [
    'gnomad_exomes_mean', 'gnomad_exomes_median', 'gnomad_exomes_over_20', 'gnomad_exomes_over_30',
    'gnomad_genomes_mean', 'gnomad_genomes_median', 'gnomad_genomes_over_20', 'gnomad_genomes_over_30',
]

# OR if using struct format:
COVERAGE_COLUMNS = ['coverage']  # Struct with exome/genome nested
```

### Phase 3: Backend Updates

**File**: `browser/backend/app.py`

Handle coverage tracks in `get_track_data()`:

```python
# In get_track_data(), add handling for coverage tracks:
elif track_id.startswith('gnomad_coverage_'):
    metric = track_id.replace('gnomad_coverage_', '')  # e.g., "over_20", "mean"
    exome_col = f'gnomad_exomes_{metric}'
    genome_col = f'gnomad_genomes_{metric}'

    values = []
    for _, row in subset.iterrows():
        values.append({
            "filtered_idx": int(row['filtered_idx']),
            "exome": row.get(exome_col) if pd.notna(row.get(exome_col)) else None,
            "genome": row.get(genome_col) if pd.notna(row.get(genome_col)) else None,
        })
    return {"track_id": track_id, "values": values}
```

### Phase 4: Track Tree Update

**File**: `browser/backend/track_tree.py`

Add coverage category with multiple overlayed tracks:

```python
TRACK_TREE = {
    "Coverage": {
        "gnomAD (Exome + Genome Overlayed)": [
            "gnomad_coverage_mean",
            "gnomad_coverage_median",
            "gnomad_coverage_over_20",
            "gnomad_coverage_over_30",
        ]
    },
    # ... existing categories
}
```

Each track (e.g., `gnomad_coverage_over_20`) renders both exome and genome data overlayed.

### Phase 5: Frontend Visualization

**File**: `browser/frontend/index.html`

#### A. Track Type Detection
```javascript
function isCoverageTrack(trackId) {
    return trackId.startsWith('gnomad_coverage_');
}
```

#### B. New Rendering Function
```javascript
function drawCoverageStacked(ctx, trackData, xScale, yScale, height, width) {
    // gnomAD-style overlayed area charts
    const exomeColor = 'rgba(70, 130, 180, 0.7)';   // Steel blue, 70% opacity
    const genomeColor = 'rgba(115, 171, 61, 0.5)';  // Green, 50% opacity

    // Draw genome first (background)
    ctx.fillStyle = genomeColor;
    ctx.beginPath();
    ctx.moveTo(0, height);
    trackData.forEach((d, i) => {
        if (d.genome && d.genome.over_20 != null) {
            const x = xScale(d.filtered_idx);
            const y = yScale(d.genome.over_20);
            if (i === 0) ctx.lineTo(x, y);
            else ctx.lineTo(x, y);
        }
    });
    ctx.lineTo(width, height);
    ctx.closePath();
    ctx.fill();

    // Draw exome on top (foreground)
    ctx.fillStyle = exomeColor;
    ctx.beginPath();
    ctx.moveTo(0, height);
    trackData.forEach((d, i) => {
        if (d.exome && d.exome.over_20 != null) {
            const x = xScale(d.filtered_idx);
            const y = yScale(d.exome.over_20);
            if (i === 0) ctx.lineTo(x, y);
            else ctx.lineTo(x, y);
        }
    });
    ctx.lineTo(width, height);
    ctx.closePath();
    ctx.fill();
}
```

#### C. Legend for Coverage Track
```javascript
function drawCoverageLegend(ctx, x, y) {
    // Exome legend
    ctx.fillStyle = 'rgba(70, 130, 180, 0.7)';
    ctx.fillRect(x, y, 15, 10);
    ctx.fillStyle = '#333';
    ctx.fillText('Exome', x + 20, y + 9);

    // Genome legend
    ctx.fillStyle = 'rgba(115, 171, 61, 0.5)';
    ctx.fillRect(x + 80, y, 15, 10);
    ctx.fillStyle = '#333';
    ctx.fillText('Genome', x + 100, y + 9);
}
```

#### D. Y-axis Scale
```javascript
// For over_X metrics: 0 to 1 (fraction of samples)
// For mean/median: 0 to max value (e.g., 100x)
const yScale = d3.scaleLinear()
    .domain([0, 1])  // For over_20 metric
    .range([height, 0]);
```

---

## Files to Modify

| Phase | File | Changes |
|-------|------|---------|
| 1 | `notebooks/browser_data.ipynb` | Add coverage join/annotation |
| 2 | `scripts/preprocess_browser_data.py` | Include coverage columns |
| 3 | `browser/backend/app.py` | Add coverage endpoint or track handler |
| 4 | `browser/backend/track_tree.py` | Add Coverage category |
| 5 | `browser/frontend/index.html` | Add `drawCoverageStacked()`, legend, track detection |

---

## Design Decisions (Confirmed)

1. **Multiple tracks per metric**: Separate tracks for mean, median, over_20, over_30
2. **Overlayed exome+genome**: Each track shows both with semi-transparent area charts
3. **gnomAD colors**: Exome = steel blue (0.7 opacity), Genome = green (0.5 opacity)
4. **Y-axis**: 0-1 for over_X metrics, auto-scale for mean/median
