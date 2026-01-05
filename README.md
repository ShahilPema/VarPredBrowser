# VarPredBrowser

A genomic variant browser with data processing pipelines, interactive genome visualization, and 3D protein structure viewing.

## Overview

This repository contains:

1. **Data Processing Pipelines** - Hail-based pipelines to merge genomic data sources (ClinVar, constraint metrics, pathogenicity scores, domain annotations) into browser-ready Parquet files.

2. **Interactive Genome Browser** - A FastAPI backend + D3.js frontend for exploring genomic data in compressed coordinate space.

3. **3D Protein Viewer** - Mol* integration for visualizing protein structures with residue-level constraint coloring.

4. **PLM vs O/E Analysis** - Tools for analyzing the relationship between predicted pathogenicity and population constraint.

## Features

- **Compressed Coordinate System**: View all coding positions across the genome with seamless navigation
- **Multi-Track Visualization**: Display constraint metrics, ClinVar annotations, pathogenicity scores, and more
- **Interactive Search**: Search by gene name, genomic position, or amino acid range
- **3D Structure Integration**: View AlphaFold structures colored by constraint/pathogenicity scores
- **Dual Data Source Support**: Works with RGC and gnomAD v4 variant data
- **Multiple Prediction Models**: Compare Core, Complete, and Constraint pathogenicity models

## Repository Structure

```
VarPredBrowser/
├── README.md                   # This file
├── requirements.txt            # Full dependencies (Hail + browser)
├── requirements-browser.txt    # Browser-only dependencies
├── .gitignore
│
├── config/
│   └── paths.yaml              # Data file paths and browser configuration
│
├── browser/                    # Interactive genome browser
│   ├── backend/
│   │   ├── app.py              # FastAPI application
│   │   ├── config.py           # Configuration loader
│   │   ├── models.py           # Pydantic models
│   │   ├── coordinate_mapper.py # Genomic-protein coordinate mapping
│   │   └── track_tree.py       # Track hierarchy definitions
│   │
│   └── frontend/
│       ├── index.html          # Main browser interface
│       ├── d3.v7.min.js        # D3.js library
│       ├── package.json        # Node.js dependencies
│       ├── vite.config.mjs     # Vite build configuration
│       └── src/
│           └── molstar-plugin.ts  # Mol* 3D viewer integration
│
├── notebooks/
│   ├── browser_data.ipynb      # Browser data preparation pipeline
│   └── plm_oe_curves.ipynb     # PLM vs O/E curve analysis
│
├── scripts/
│   ├── preprocess_browser_data.py  # Generate axis tables from merged data
│   ├── build_protein_map.py        # Build genomic-protein coordinate maps
│   └── fetch_alphafold.py          # Fetch AlphaFold structures
│
├── docs/
│   ├── INPUT_FILES.md          # Input data documentation
│   └── METHODS.md              # Methodology documentation
│
├── data/                       # Processed browser data (gitignored)
│   ├── any_count_gt0.parquet   # All coding positions
│   ├── mis_count_gt0.parquet   # Missense-possible positions
│   ├── gene_index_*.parquet    # Gene lookup tables
│   └── structures/             # Protein structure files
│
└── output/                     # Analysis outputs (gitignored)
    ├── figures/
    └── curve_data/
```

## Quick Start

### Option 1: Browser Only (Pre-processed Data)

If you have pre-processed data, you only need to run the browser:

```bash
# Install browser dependencies
pip install -r requirements-browser.txt

# Install frontend dependencies
cd browser/frontend
npm install
npm run build
cd ../..

# Start the server
python -m browser.backend.app
```

Access the browser at: http://localhost:8000

### Option 2: Full Pipeline

For processing raw data and running the browser:

```bash
# Install all dependencies
pip install -r requirements.txt

# Configure data paths
# Edit config/paths.yaml with your data file locations

# Step 1: Generate merged browser data (requires Hail)
# Run notebooks/browser_data.ipynb

# Step 2: Preprocess for browser
python scripts/preprocess_browser_data.py

# Step 3 (Optional): Build protein maps for specific genes
python scripts/build_protein_map.py --gene SCN2A

# Step 4 (Optional): Fetch AlphaFold structures
python scripts/fetch_alphafold.py --gene SCN2A

# Step 5: Start the browser
python -m browser.backend.app
```

## Data Pipeline

```
[Hail Pipeline]                    [Preprocessing]                  [Browser]
notebooks/browser_data.ipynb  -->  preprocess_browser_data.py  --> browser/backend/app.py
       |                                    |                              |
       v                                    v                              v
rgc_browser_data_merged.parquet      data/*.parquet              FastAPI + D3.js
(33.3M rows, 278 cols)               (axis tables)               http://localhost:8000
```

## Configuration

Edit `config/paths.yaml` to configure:

```yaml
data_sources:
  browser_data: "/path/to/rgc_browser_data_merged.parquet"
  # ... other Hail table paths

browser:
  data_dir: "data"
  port: 8000
```

Environment variables can override config:
- `VARPRED_DATA_DIR`: Browser data directory
- `VARPRED_PORT`: Server port

## API Endpoints

The browser exposes a RESTful API:

| Endpoint | Description |
|----------|-------------|
| `GET /api/health` | Health check |
| `GET /api/filters` | List available data filters |
| `GET /api/track-tree` | Get track hierarchy |
| `GET /api/filtered-window` | Get positions in coordinate window |
| `GET /api/track-data` | Get track values for window |
| `GET /api/search/gene` | Search by gene symbol |
| `GET /api/search/gene/autocomplete` | Gene autocomplete |
| `GET /api/protein/{gene}` | Get protein info |
| `GET /api/protein/{gene}/residue-scores` | Get per-residue scores |
| `GET /api/structure/{gene}` | Get structure metadata |
| `GET /api/structure/{gene}/file/{type}` | Download structure file |

Full API documentation available at: http://localhost:8000/docs

## Development

### Frontend Development (with hot reload)

```bash
# Terminal 1: Start backend
python -m browser.backend.app

# Terminal 2: Start frontend dev server
cd browser/frontend
npm run dev
```

Access at: http://localhost:5173 (proxies API to backend)

### Production Build

```bash
cd browser/frontend
npm run build  # Outputs to dist/

# Server will automatically serve from dist/
python -m browser.backend.app
```

## Data Requirements

For the full pipeline, you need access to:

| File | Size | Description |
|------|------|-------------|
| `rgc_browser_data_merged.parquet` | ~21GB | Pre-merged browser data |
| `rgc_scaled.ht` | ~100GB | RGC variant data (optional) |
| `AOU_RGC_All_preds.ht` | ~50GB | ML predictions (optional) |

See [docs/INPUT_FILES.md](docs/INPUT_FILES.md) for detailed documentation.

## Prediction Models

Three pathogenicity prediction models are visualized:

- **Constraint**: Uses population constraint features
- **Core**: Sequence-based features only
- **Complete**: Full feature set including AlphaFold structural features

## License

[License information to be added]
