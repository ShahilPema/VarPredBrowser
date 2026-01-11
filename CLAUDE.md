# CLAUDE.md - AI Assistant Guide for VarPredBrowser

## Project Overview

VarPredBrowser is a genomic variant browser with integrated 3D protein structure visualization. It combines:

- **Hail-based data pipelines** - Merge genomic data (ClinVar, constraint metrics, pathogenicity scores, domain annotations)
- **FastAPI backend** - Serve processed data via REST API with coordinate mapping
- **D3.js frontend** - Interactive genome visualization in compressed coordinate space
- **Mol* 3D viewer** - Protein structure visualization with constraint/pathogenicity coloring

## Architecture

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                              DATA PIPELINE                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│  notebooks/browser_data.ipynb                                                 │
│    ↓ (Hail processing)                                                       │
│  rgc_browser_data_merged.parquet (33M rows × 278 cols)                       │
│    ↓                                                                         │
│  scripts/preprocess_browser_data.py                                          │
│    ↓                                                                         │
│  data/*.parquet (axis tables + gene indexes)                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                               BROWSER                                         │
├──────────────────────────────────────────────────────────────────────────────┤
│  browser/backend/                     browser/frontend/                       │
│  ├── app.py (FastAPI)                 ├── index.html (D3.js visualization)   │
│  ├── config.py                        ├── src/molstar-plugin.ts (3D viewer)  │
│  ├── models.py (Pydantic)             ├── package.json (Node deps)           │
│  ├── coordinate_mapper.py             └── vite.config.mjs                    │
│  └── track_tree.py                                                           │
└──────────────────────────────────────────────────────────────────────────────┘
```

## Key Files

### Backend (`browser/backend/`)

| File | Purpose |
|------|---------|
| `app.py` | FastAPI application with all API endpoints. Loads axis tables into memory on startup. |
| `config.py` | Configuration loader from `config/paths.yaml` with env var overrides |
| `models.py` | Pydantic data models for API request/response |
| `coordinate_mapper.py` | Genomic-to-protein coordinate conversion, variant data extraction |
| `track_tree.py` | Hierarchical track definitions (O/E ratios, VIRs, ClinVar, predictions) |

### Frontend (`browser/frontend/`)

| File | Purpose |
|------|---------|
| `index.html` | Main browser interface (~298KB) with D3.js visualization code |
| `src/molstar-plugin.ts` | Mol* 3D protein viewer integration with sphere halo highlighting |
| `package.json` | Node.js deps: molstar ^5.5.0, react ^19.2.3, vite ^7.3.0 |
| `vite.config.mjs` | Build config, proxies to backend on port 8000 |

### Data Pipeline

| File | Purpose |
|------|---------|
| `notebooks/browser_data.ipynb` | Main Hail pipeline - merges all data sources |
| `scripts/preprocess_browser_data.py` | Generates axis tables from merged parquet |
| `scripts/build_protein_map.py` | Creates genomic-to-protein coordinate maps |
| `scripts/fetch_alphafold.py` | Downloads AlphaFold structures |
| `scripts/fetch_interpro_domains.py` | Fetches domain annotations via InterPro API |

### Configuration

| File | Purpose |
|------|---------|
| `config/paths.yaml` | Data source paths, Hail settings, browser config |
| `requirements.txt` | Full deps (Hail + browser) |
| `requirements-browser.txt` | Browser-only deps (no Hail) |

## Development Commands

### Run Backend
```bash
python -m browser.backend.app
# Runs on http://localhost:8000
```

### Run Frontend (Development with hot reload)
```bash
cd browser/frontend
npm install  # First time only
npm run dev
# Runs on http://localhost:5173, proxies API to :8000
```

### Build Frontend for Production
```bash
cd browser/frontend
npm run build  # Outputs to dist/
```

### Run Data Preprocessing
```bash
python scripts/preprocess_browser_data.py --input /path/to/merged.parquet --output data/
```

## Planned Updates

Seven plan files in `plans/` describe upcoming features:

### 1. browser-data-refactor.md (Priority: High)
- Convert tuple aggregations to structs for better Parquet/JSON compatibility
- Affects: `preds.Constraint`, `dbnsfp_stacked.AlphaMissense`, `variant_consequences`
- Downstream: Update `coordinate_mapper.py` to use `.alt`, `.pred`, `.score` instead of `[0]`, `[1]`

### 2. variant-protein-overlay.md (Priority: High)
- Overlay variant tracks (ClinVar, Training, Constraint) as spheres on 3D structure
- Max 2 simultaneous overlays: blue (#0066CC), cyan (#00CED1)
- UI: Button below flip button on track colorbar

### 3. add-gnomad-constraint.md (Priority: Medium)
- Add gnomAD O/E constraint metrics parallel to RGC
- Source: `constraint_metrics_by_locus_revft.ht`
- Naming: `gnomad_mis_21bp_oe` (simplified vs RGC)

### 4. coverage-track.md (Priority: Medium)
- gnomAD exome/genome coverage visualization
- Area charts (large regions), bar charts (<100bp)
- Data already in: `gnomad_exomes_over_20`, `gnomad_genomes_over_20`

### 5. variant-frequency-track.md (Priority: Medium)
- Log-scaled ellipse encoding for allele frequency
- Cohort colors: blue (exomes), green (genomes), purple (joint), orange (RGC)

### 6. clinvar-track-upgrade.md (Priority: Low)
- Modernize ClinVar track with new data structure
- Add star filtering, consequence filtering, expand/collapse views

### 7. miscellaneous.md (Priority: Mixed)
- Yellow sphere highlight fix (click outside to deselect)
- RGC tracks hierarchy flattening (5 levels → 3)
- 3D panel show/hide tab
- Track reordering via drag-and-drop
- Remove title bar for more space

## Data Structures

### Track Types

```python
# From track_tree.py
CONSTRAINT_STACKED_FIELDS = {'Constraint', 'Core', 'Complete'}
DBNSFP_STACKED_FIELDS = {'AlphaMissense_stacked', 'ESM1b_stacked'}
```

### Stacked Track Data (Current - Tuples)
```python
# Constraint predictions: array of (allele, pred, n_pred)
row['Constraint'] = [('T', 0.85, 3), ('A', 0.72, 3)]

# dbNSFP stacked: array of (allele, score, percentile)
row['AlphaMissense_stacked'] = [('T', 0.95, 98.2), ('A', 0.45, 52.1)]
```

### Stacked Track Data (Planned - Structs)
```python
# After refactoring
row['Constraint'] = [{'alt': 'T', 'pred': 0.85, 'n_pred': 3}, ...]
row['AlphaMissense_stacked'] = [{'alt': 'T', 'score': 0.95, 'percentile': 98.2}, ...]
```

### Filter Definitions
```python
FILTERS = {
    'mis_count_gt0': 'Positions where missense variants are possible',
    'any_count_gt0': 'All coding positions where any variant is possible'
}
```

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /api/health` | Health check, shows loaded filters |
| `GET /api/filters` | Available data filters |
| `GET /api/track-tree` | Hierarchical track tree for UI |
| `GET /api/filtered-window` | Positions in compressed coordinate window |
| `GET /api/track-data` | Track values for coordinate window |
| `GET /api/search/gene` | Gene symbol lookup |
| `GET /api/search/gene/autocomplete` | Gene autocomplete |
| `GET /api/protein/{gene}` | Protein info (length, domains, structure) |
| `GET /api/protein/{gene}/residue-scores` | Per-residue aggregated scores |
| `GET /api/structure/{gene}` | Structure metadata |
| `GET /api/structure/{gene}/file/{type}` | Structure file (PDB) download |

## Coding Patterns

### Backend Data Extraction
```python
# From coordinate_mapper.py
def extract_constraint_variants(raw_value):
    """Extract (allele, pred) tuples from constraint stacked field."""
    # Returns list of (allele, prediction) for tooltip display

def extract_dbnsfp_stacked_variants(raw_value):
    """Extract (allele, score, percentile) from dbNSFP stacked field."""
```

### Frontend Track Classification
```javascript
// From index.html
function isClinVarStackedTrack(trackId) { return trackId === 'clinvar.clinvar_label_list'; }
function isTrainingTrack(trackId) { return trackId.startsWith('training.train_counts.'); }
function isConstraintStackedTrack(trackId) { return ['Constraint', 'Core', 'Complete'].includes(trackId); }
```

### Mol* Integration
```typescript
// From molstar-plugin.ts
export async function createSphereHalo(loci, residueNum, sizeFactor, alpha) {
    // Creates yellow sphere highlight around clicked residue
}
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `VARPRED_DATA_DIR` | `data/` | Browser data directory |
| `VARPRED_BIGWIG_DIR` | (from config) | BigWig files directory |
| `VARPRED_HOST` | `0.0.0.0` | Server host |
| `VARPRED_PORT` | `8000` | Server port |

## Key Dependencies

### Python (Backend)
- `fastapi>=0.104.1` - Web framework
- `polars>=0.20.0` - DataFrame operations (faster than pandas)
- `pyarrow>=14.0.1` - Parquet I/O
- `pybigwig>=0.3.22` - BigWig track files
- `hail>=0.2.134` - Genomic data processing (pipeline only)

### Node.js (Frontend)
- `molstar@^5.5.0` - 3D protein structure viewer
- `react@^19.2.3` - UI components for Mol*
- `vite@^7.3.0` - Build tool with HMR
- `typescript@^5.9.3` - Type checking

## Common Tasks

### Adding a New Track
1. Add field to `browser_data.ipynb` aggregation
2. Re-run `preprocess_browser_data.py`
3. Add to `track_tree.py` hierarchy
4. If stacked track: add to `CONSTRAINT_STACKED_FIELDS` or `DBNSFP_STACKED_FIELDS`

### Adding a New API Endpoint
1. Add route in `browser/backend/app.py`
2. Add Pydantic model in `models.py` if needed
3. Test at http://localhost:8000/docs

### Modifying 3D Viewer
1. Edit `browser/frontend/src/molstar-plugin.ts`
2. Rebuild: `cd browser/frontend && npm run build`

### Fixing Coordinate Mapping
1. Check `coordinate_mapper.py` for extraction logic
2. Verify protein map exists in `data/{gene}_protein_map.parquet`

## Testing

No formal test suite exists yet. Manual testing workflow:

1. Start backend: `python -m browser.backend.app`
2. Start frontend: `cd browser/frontend && npm run dev`
3. Open http://localhost:5173
4. Test: gene search, track selection, 3D structure loading, residue clicking

## Data Sources (from config/paths.yaml)

| Source | Path | Description |
|--------|------|-------------|
| `rgc_scaled` | `/storage/zoghbi/.../rgc_scaled.ht` | RGC variant data (817M variants) |
| `gnomad_scaled` | `/storage/zoghbi/.../gnomadV4_scaled.ht` | gnomAD v4 variant data |
| `predictions` | `/local/.../AOU_RGC_All_preds.ht` | ML predictions (Constraint, Core, Complete) |
| `browser_data` | `/storage/zoghbi/.../rgc_browser_data_merged.parquet` | Pre-merged browser data |
| `dbnsfp` | `/local/.../dbnsfp/...` | dbNSFP annotation scores |
| `gnomad_exomes_coverage` | `/storage/zoghbi/.../gnomadV4_exomes_coverage_struct.ht` | Coverage data |

## Notes for AI Assistants

1. **Frontend is large**: `index.html` is ~298KB with embedded D3 code. Consider reading in chunks if needed.

2. **Tuple → Struct migration**: Current stacked data uses tuples (`[0]`, `[1]` indexing). Plans call for structs (`.alt`, `.pred`, `.score`). Check which format is active before modifying extraction code.

3. **Filter reactivity**: When modifying track filters (ClinVar labels, etc.), consider if 3D overlay spheres need updating.

4. **Coordinate systems**: Three coordinate systems exist:
   - Genomic (chr:pos)
   - Compressed/filtered (filtered_idx)
   - Protein (residue number)

5. **Configuration priority**: Environment variables override `config/paths.yaml` values.

6. **Hail vs Polars**: Pipeline uses Hail (Spark-based), browser uses Polars (much faster for single-machine).

7. **Hail table schema inspection**: Every Hail table directory contains a `metadata.json.gz` file at its root. Use this to inspect the schema without loading the table into Hail:
   ```bash
   zcat /path/to/table.ht/metadata.json.gz | python3 -m json.tool
   ```
   This is especially useful for discovering field names and types when integrating new data sources. For example, the gnomAD constraint metrics use field pattern `gnomad_{consequence}_{cohort}_{window}_{metric}` (e.g., `gnomad_mis_exomes_21bp_oe`).

8. **Hail key field modification**: You cannot overwrite key fields using `annotate`, `select`, or `drop`. If you get the error `ExpressionException: 'Table.select': cannot overwrite key field 'locus' with annotate, select or drop; use key_by to modify keys`:
   ```python
   # ALWAYS check key fields first before using select()
   key_fields = list(ht.key)
   print(f"Key fields (auto-included, cannot be in select): {key_fields}")

   # Wrong - will error if 'locus' or 'region' are key fields
   ht = ht.select('locus', 'region', 'other_field')

   # Correct - exclude key fields from select list, they're auto-included
   cols_to_keep = [col for col in my_cols if col not in key_fields]
   ht = ht.select(*cols_to_keep)

   # If you need to modify a key field, unkey first
   ht = ht.key_by()  # Remove all keys
   ht = ht.annotate(locus = new_locus)
   ht = ht.key_by('locus', 'alleles')  # Re-key
   ```

9. **Hail join optimization**: When merging datasets, join once to get the struct, then annotate fields from it. Each table lookup is a join operation, so multiple lookups to the same table means multiple joins:
   ```python
   # Wrong - 4 separate join operations
   ht = ht.annotate(
       field1 = other_ht[ht.key].field1,
       field2 = other_ht[ht.key].field2,
       field3 = other_ht[ht.key].field3,
       field4 = other_ht[ht.key].field4,
   )

   # Correct - 1 join operation, then extract fields
   ht = ht.annotate(_other = other_ht[ht.key])
   ht = ht.annotate(
       field1 = ht._other.field1,
       field2 = ht._other.field2,
       field3 = ht._other.field3,
       field4 = ht._other.field4,
   )
   ht = ht.drop('_other')
   ```

10. **Always check table keys before joining**: Tables may have empty keys or unexpected key structure. A join like `ht[other_ht.locus]` requires `ht` to be keyed by `locus`:
    ```python
    # ALWAYS check keys before joining
    print(f"Table key: {list(ht.key) or '<<<EMPTY KEY>>>'}")

    # If empty or wrong key, must key_by first
    ht = ht.key_by('locus')  # Now can join by locus
    ```
    Common error: `ExpressionException: Key type mismatch: Table key: <<<empty key>>>`

11. **Testing Hail code**: Before modifying notebook cells, create a test script to validate the logic:
    ```python
    # scripts/test_hail_logic.py
    import hail as hl
    hl.init(master='local[4]', quiet=True)

    # Load table and inspect key/row structure
    ht = hl.read_table('/path/to/table.ht')
    print(f"Key: {list(ht.key)}")
    print(f"Row fields: {list(ht.row)[:10]}...")

    # Test operations before putting in notebook
    ```
    Run with: `python scripts/test_hail_logic.py 2>&1`

12. **Locus-based vs variant-based tables**: The base browser table is keyed by `locus` (position-level), while variant data like allele frequencies are keyed by `(locus, alleles)`. You cannot directly join these - you must aggregate variant-level data by locus first:
    ```python
    # Wrong - base_ht has no 'alleles' field
    base_ht = base_ht.annotate(
        af = variant_ht[base_ht.locus, base_ht.alleles].af  # Error!
    )

    # Correct - aggregate variant data by locus into array of structs
    variant_by_locus = variant_ht.group_by('locus').aggregate(
        variants = hl.agg.collect(hl.struct(
            alt = variant_ht.alleles[1],
            af = variant_ht.info.AF,
            ac = variant_ht.info.AC,
        ))
    )
    base_ht = base_ht.annotate(
        variants = variant_by_locus[base_ht.locus].variants
    )
    ```
    This pattern is used for allele frequencies (RGC, gnomAD) and matches the stacked score format (`AlphaMissense_stacked`, `ESM1b_stacked`).
