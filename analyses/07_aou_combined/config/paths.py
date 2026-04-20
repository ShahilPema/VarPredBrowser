"""Single-edit-point path config.

Defaults derive ROOT from this file's own location, so the dir works
unchanged on BCM, the AoU workbench, or anywhere else cloned. Override
via VARPRED_07_ROOT / VARPRED_07_INPUTS env vars only if you want to
relocate inputs/outputs out of the analysis dir.
"""
import os
from pathlib import Path

# -----------------------------------------------------------------------------
# Local roots — derived from this file's location (config/paths.py is at
# <ROOT>/config/paths.py, so parent.parent is the analysis dir).
# -----------------------------------------------------------------------------

_HERE = Path(__file__).resolve().parent.parent

ROOT = Path(os.environ.get('VARPRED_07_ROOT', str(_HERE)))
INPUTS_DIR = Path(os.environ.get('VARPRED_07_INPUTS', str(ROOT / 'inputs')))
OUTPUT_DIR = ROOT / 'output'
PER_SITE_DIR = OUTPUT_DIR / 'per_site'
CURVES_DIR = OUTPUT_DIR / 'curves'
FIG_DIR = OUTPUT_DIR / 'figures'
LOGS_DIR = OUTPUT_DIR / 'logs'

for d in (INPUTS_DIR, OUTPUT_DIR, PER_SITE_DIR, CURVES_DIR, FIG_DIR, LOGS_DIR):
    d.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# GCS
# -----------------------------------------------------------------------------

GCS_BUCKET = 'gs://zoghbi-lab-data/shahil'
GCS_INPUTS_PREFIX = f'{GCS_BUCKET}/07_aou_combined_inputs'
GCS_OUTPUTS_PREFIX = f'{GCS_BUCKET}/07_aou_combined_outputs'

# -----------------------------------------------------------------------------
# Shipped input files (consumed by 01_build_per_site.py on workbench)
# -----------------------------------------------------------------------------

BASE_TABLE_PARQUET = INPUTS_DIR / 'base_table.parquet'
INPUT_MANIFEST = INPUTS_DIR / 'MANIFEST.json'

# -----------------------------------------------------------------------------
# Per-site outputs (one per dataset)
# -----------------------------------------------------------------------------

DATASETS = ('gnomad_only', 'aou_only', 'combined')

def per_site_path(dataset: str) -> Path:
    assert dataset in DATASETS, dataset
    return PER_SITE_DIR / f'per_site_{dataset}.parquet'


# -----------------------------------------------------------------------------
# Score tag definitions (pared down per cohort comparison scope)
# -----------------------------------------------------------------------------

# `raw`            : in-house Full_Raw_Core (best in-house performer in last eval)
# `popeve_neg`     : external popEVE
# `alphamissense`  : external AlphaMissense (top external in last eval)
# Quadprog runs on in-house tags only (see TAGS_PER_METHOD in _fit_common.py),
# so pruning to one in-house tag means quadprog now runs on just `raw`.
INHOUSE_TAGS = ('raw',)
EXTERNAL_TAGS = ('popeve_neg', 'alphamissense')
ALL_TAGS = INHOUSE_TAGS + EXTERNAL_TAGS

TAG_TO_SCORE_COL = {
    'raw':           'score_raw',
    'popeve_neg':    'score_popeve_neg',
    'alphamissense': 'score_alphamissense',
}

# -----------------------------------------------------------------------------
# AoU CDR pin (override per release)
# -----------------------------------------------------------------------------

AOU_CDR_VERSION = os.environ.get('AOU_CDR_VERSION', 'v8')
