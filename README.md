# VarPredBrowser

A genomic variant browser data pipeline and PLM (Protein Language Model) vs O/E (Observed/Expected) curve analysis toolkit.

## Overview

This repository contains tools for:

1. **Browser Data Preparation** - Merge multiple genomic data sources (ClinVar, constraint metrics, pathogenicity scores, domain annotations) into browser-ready Parquet files for the Gosling genome visualization platform.

2. **PLM vs O/E Curve Analysis** - Generate per-gene curves showing the relationship between predicted pathogenicity (PLM scores) and population constraint (O/E ratios).

## Features

- **Dual Data Source Support**: Generate curves from either RGC (Regeneron Genetics Center) or gnomAD v4 variant data
- **Multiple Prediction Models**: Compare Core, Complete, and Constraint pathogenicity models
- **Population Constraint Metrics**: Sliding window O/E ratios ordered by prediction score
- **ClinVar Validation**: Overlay pathogenic variants on curves
- **Configurable Paths**: Easy configuration via YAML file

## Repository Structure

```
VarPredBrowser/
├── README.md                   # This file
├── requirements.txt            # Python dependencies
├── .gitignore                  # Git ignore patterns
│
├── config/
│   └── paths.yaml              # Data file paths configuration
│
├── notebooks/
│   ├── browser_data.ipynb      # Browser data preparation pipeline
│   └── plm_oe_curves.ipynb     # PLM vs O/E curve analysis
│
├── scripts/
│   └── preprocess_browser_data.py  # Browser data preprocessing
│
├── docs/
│   ├── INPUT_FILES.md          # Input data documentation
│   └── METHODS.md              # Methodology documentation
│
└── output/                     # Generated outputs (gitignored)
    ├── figures/
    └── curve_data/
```

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Configure Data Paths

Edit `config/paths.yaml` to point to your data files:

```yaml
data_sources:
  rgc_scaled: "/path/to/rgc_scaled.ht"
  gnomad_scaled: "/path/to/gnomadV4_scaled.ht"
  predictions: "/path/to/AOU_RGC_All_preds.ht"
  clinvar: "/path/to/clinvar.ht"
```

### 3. Run PLM vs O/E Analysis

Open `notebooks/plm_oe_curves.ipynb` and configure:

```python
# Select data source: 'rgc' or 'gnomad'
DATA_SOURCE = 'rgc'

# Select prediction model: 'constraint', 'core', or 'complete'
PREDICTION_MODEL = 'constraint'
```

Then run all cells to generate per-gene curves.

## Data Requirements

This pipeline requires access to large Hail tables stored externally:

| File | Size | Description |
|------|------|-------------|
| `rgc_scaled.ht` | ~100GB | RGC variant data with mutation rates |
| `gnomadV4_scaled.ht` | ~150GB | gnomAD v4 variant data |
| `AOU_RGC_All_preds.ht` | ~50GB | ML predictions (Core/Complete/Constraint) |
| `cv_38_final.ht` | ~2GB | ClinVar pathogenic variants |

See [docs/INPUT_FILES.md](docs/INPUT_FILES.md) for detailed data documentation.

## Prediction Models

Three pathogenicity prediction models are available:

- **Constraint**: Uses population constraint features
- **Core**: Sequence-based features only (no structural data)
- **Complete**: Full feature set including AlphaFold structural features

## Key Outputs

### PLM vs O/E Curves

For each gene, the notebook generates:
- Individual gene curves (PLM score on x-axis, local O/E on y-axis)
- ClinVar pathogenic variant overlay
- Curve category classification (monotonic decreasing, flat, other)
- Gene-level statistics (correlation, O/E range)

### Browser Data

The browser data pipeline produces:
- Locus-level Parquet files for Gosling visualization
- Gene index files for navigation
- ~33 million positions with ~280 annotation columns

## Methodology

See [docs/METHODS.md](docs/METHODS.md) for details on:
- Sliding window O/E calculation
- PLM score ordering approach
- Edge handling and window expansion
- Statistical validation

## Citation

If you use this code, please cite:

[Citation information to be added]

## License

[License information to be added]
