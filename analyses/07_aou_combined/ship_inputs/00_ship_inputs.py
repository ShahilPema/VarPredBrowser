#!/usr/bin/env python3
"""Build a single consolidated base_table.parquet at BCM and push to GCS.

This is the ONLY script in 07_aou_combined that runs at BCM. Everything else
runs on the AoU workbench and consumes this one parquet.

Output schema (one row per (locus, alleles, transcript) for missense+synonymous
sites with a defined roulette mutation rate):

  KEYS         : locus_str, alleles_str, transcript
  metadata     : gene_symbol, most_deleterious_consequence_cds,
                 is_chrX_nonPAR, roulette_AR_MR, roulette_syn_scaling_site,
                 PAR_X
  gnomAD obs   : gnomad_in_table, gnomad_ac, gnomad_is_pass
  filters      : f_lowcov, f_ab_low                 (gnomAD-derived)
  in-house sc  : score_c_iso, score_raw, score_c_noqvd_ps, score_core_ps
  external sc  : score_popeve_neg, score_eve_neg, score_esm1v_neg,
                 score_alphamissense

Drops vs the BCM `base_all_joins.ht`:
  RGC fields, plm_score probe, somatic flags (f_somatic, f_ab2qd_*),
  gnomAD/RGC overfilter z-score / pass93 / gene_no / synoe flags, all rows
  whose consequence is not missense or synonymous.

Then writes MANIFEST.json with row count + per-partition md5, and pushes both
to gs://zoghbi-lab-data/shahil/07_aou_combined_inputs/.

Run at BCM with: python3 -u ship_inputs/00_ship_inputs.py
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

# -----------------------------------------------------------------------------
# Paths (BCM-side; do not import config/paths.py which is workbench-targeted)
# -----------------------------------------------------------------------------

ALL_SITES_TABLE = '/storage/zoghbi/home/u235147/merged_vars/all_sites_annotated_4.ht'
GNOMAD_SITES_TABLE = '/local/tmp/gnomad_qc_checkpoints/gnomad.exomes.v4.1.sites.ht'
SCORE_INHOUSE_TABLE = '/local/blake_checkpoints/checkpoints/Inference/Predictions/AOU_RGC_with_perc_and_probs.ht'
ENST_HGNC_TABLE = '/storage/zoghbi/data/sharing/hail_tables/ENST_HGNC_mapping.ht'
FILTER_FLAGS_PARQUET = '/storage/zoghbi/home/u235147/merged_vars/independent_filter_analysis/results/position_filter_flags.parquet'

# In-house score columns on SCORE_INHOUSE_TABLE
INHOUSE_SCORE_FIELDS = {
    'score_c_iso':      'Full_Core_AOUN4_No_aapos_pred',
    'score_raw':        'Full_Raw_Core_AOUN4_NA_ps_pred',
    'score_c_noqvd_ps': 'Full_Core_no_qvd_AOUN4_NA_ps_pred',
    'score_core_ps':    'Full_Core_AOUN4_NA_ps_pred',
}

# External scores live in scores_extra.parquet; reuse it as-is.
SCORES_EXTRA_PARQUET = '/storage/zoghbi/home/u235147/VarPredBrowser/analyses/04_bernoulli_fitter/output/scores_extra.parquet'
EXTERNAL_SCORE_COLS = (
    'score_popeve_neg', 'score_eve_neg', 'score_esm1v_neg', 'score_alphamissense',
)

ROOT = Path(__file__).resolve().parent.parent
STAGE_DIR = ROOT / 'inputs_stage'
STAGE_DIR.mkdir(parents=True, exist_ok=True)
BASE_OUT = STAGE_DIR / 'base_table.parquet'
MANIFEST_OUT = STAGE_DIR / 'MANIFEST.json'

GCS_DEST = 'gs://zoghbi-lab-data/shahil/07_aou_combined_inputs'

SPARK_TMP = '/local/tmp/aou_combined_ship'


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# -----------------------------------------------------------------------------
# Hail base table build
# -----------------------------------------------------------------------------

def init_hail():
    import hail as hl
    cpus = int(os.environ.get('SHIP_CPUS', '110'))
    memory = int(3600 * cpus / 256)
    hl.init(
        spark_conf={
            'spark.driver.memory': f'{memory}g',
            'spark.executor.memory': f'{memory}g',
            'spark.local.dir': SPARK_TMP,
            'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={SPARK_TMP}',
            'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={SPARK_TMP}',
            'spark.ui.enabled': 'false',
        },
        master=f'local[{cpus}]', tmp_dir=SPARK_TMP, local_tmpdir=SPARK_TMP,
        quiet=True,
    )
    hl.default_reference('GRCh38')
    return hl


def build_base_table_hail():
    os.makedirs(SPARK_TMP, exist_ok=True)
    hl = init_hail()

    log('Reading source tables...')
    all_sites = hl.read_table(ALL_SITES_TABLE)
    all_sites = all_sites.filter(
        hl.is_defined(all_sites.roulette_AR_MR) & (all_sites.roulette_AR_MR > 0)
    )

    # Slim to missense + synonymous only — saves ~3-4x on output size.
    all_sites = all_sites.filter(
        (all_sites.most_deleterious_consequence_cds == 'missense_variant') |
        (all_sites.most_deleterious_consequence_cds == 'synonymous_variant')
    )
    log(f'  all_sites (missense+syn, defined MR): {all_sites.count():,}')

    gnomad = hl.read_table(GNOMAD_SITES_TABLE)
    score_inhouse = hl.read_table(SCORE_INHOUSE_TABLE).key_by('locus', 'alleles')
    enst = hl.read_table(ENST_HGNC_TABLE)

    # gnomAD per-site: pull AC + PASS, default missing -> 0/False/Not in table
    log('Joining gnomAD AC + PASS...')
    _g = gnomad[all_sites.locus, all_sites.alleles]
    ht = all_sites.annotate(
        transcript=all_sites.region.split('-')[0],
        gnomad_in_table=hl.is_defined(_g),
        gnomad_ac=hl.or_else(_g.freq[0].AC, 0),
        gnomad_is_pass=hl.or_else(hl.len(_g.filters) == 0, False),
        is_chrX_nonPAR=hl.if_else(
            (all_sites.locus.contig == 'chrX') & (hl.or_else(all_sites.PAR_X, 0) == 0),
            True, False,
        ),
    )

    log('Joining ENST -> HGNC for gene_symbol...')
    ht = ht.annotate(gene_symbol=enst[ht.transcript].HGNC)

    log('Joining in-house scores...')
    _s = score_inhouse[ht.locus, ht.alleles]
    ht = ht.annotate(**{
        out_col: _s[hail_field]
        for out_col, hail_field in INHOUSE_SCORE_FIELDS.items()
    })

    log('Joining gnomAD-derived position filter flags...')
    if os.path.exists(FILTER_FLAGS_PARQUET):
        from pyspark.sql import SparkSession
        spark_session = SparkSession.builder.getOrCreate()
        ff_df = spark_session.read.parquet(FILTER_FLAGS_PARQUET)
        ff_ht = hl.Table.from_spark(ff_df)
        ff_fields = list(ff_ht.row)
        if 'chrom' in ff_fields:
            ff_ht = ff_ht.annotate(
                locus=hl.locus(ff_ht.chrom, ff_ht.pos, reference_genome='GRCh38'))
        else:
            ff_ht = ff_ht.annotate(
                locus=hl.locus(ff_ht.contig, ff_ht.pos, reference_genome='GRCh38'))
        ff_ht = ff_ht.key_by('locus')
        # Only keep f_lowcov and f_ab_low (drop overfilter, RGC, somatic flags).
        _ff = ff_ht[ht.locus]
        ht = ht.annotate(
            f_lowcov=hl.or_else(_ff.f_lowcov, False),
            f_ab_low=hl.or_else(_ff.f_ab_low, False),
        )
    else:
        log(f'  WARNING: {FILTER_FLAGS_PARQUET} missing; setting filters to False')
        ht = ht.annotate(f_lowcov=False, f_ab_low=False)

    # Final select — only the columns the workbench will use.
    out_fields = [
        'transcript', 'gene_symbol', 'most_deleterious_consequence_cds',
        'is_chrX_nonPAR', 'PAR_X',
        'roulette_AR_MR', 'roulette_syn_scaling_site',
        'gnomad_in_table', 'gnomad_ac', 'gnomad_is_pass',
        'f_lowcov', 'f_ab_low',
    ] + list(INHOUSE_SCORE_FIELDS.keys())

    ht = ht.annotate(
        locus_str=hl.str(ht.locus),
        alleles_str=hl.str(ht.alleles),
    )
    ht = ht.key_by().select('locus_str', 'alleles_str', *out_fields)

    log(f'Base table row count: {ht.count():,}')
    log(f'Writing partitioned parquet -> {BASE_OUT}')
    ht.to_spark().write.mode('overwrite').parquet(str(BASE_OUT))

    hl.stop()
    log('Hail stage done.')


# -----------------------------------------------------------------------------
# Join external scores via Polars (already in parquet form)
# -----------------------------------------------------------------------------

def join_external_scores():
    """Left-join scores_extra into base_table via Polars (no Hail needed)."""
    import polars as pl
    log('Joining external scores via Polars...')
    base = pl.scan_parquet(f'{BASE_OUT}/*.parquet')
    extra = (pl.scan_parquet(SCORES_EXTRA_PARQUET)
               .select(['locus_str', 'alleles_str', 'transcript', *EXTERNAL_SCORE_COLS]))
    merged = base.join(extra, on=['locus_str', 'alleles_str', 'transcript'], how='left')
    tmp_out = STAGE_DIR / 'base_table_with_extra.parquet'
    merged.sink_parquet(str(tmp_out), compression='zstd', compression_level=3)
    log(f'  wrote {tmp_out}')

    # Replace base_table.parquet with the augmented one
    import shutil
    shutil.rmtree(BASE_OUT, ignore_errors=True)
    if tmp_out.is_dir():
        shutil.move(str(tmp_out), str(BASE_OUT))
    else:
        # sink_parquet wrote a single file
        BASE_OUT.mkdir(exist_ok=True)
        shutil.move(str(tmp_out), str(BASE_OUT / 'part-0.parquet'))
    log(f'  final base_table at {BASE_OUT}')


# -----------------------------------------------------------------------------
# Manifest + GCS push
# -----------------------------------------------------------------------------

def md5_file(path: Path, chunk=8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def write_manifest():
    import polars as pl
    log('Writing MANIFEST.json...')
    n_rows = pl.scan_parquet(f'{BASE_OUT}/*.parquet').select(pl.len()).collect()['len'][0]
    parts = sorted(Path(BASE_OUT).glob('*.parquet'))
    manifest = {
        'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'base_table_dir': str(BASE_OUT.name),
        'n_rows': int(n_rows),
        'n_partitions': len(parts),
        'partitions': [
            {'name': p.name, 'size_bytes': p.stat().st_size, 'md5': md5_file(p)}
            for p in parts
        ],
    }
    with open(MANIFEST_OUT, 'w') as f:
        json.dump(manifest, f, indent=2)
    log(f'  rows={n_rows:,}  partitions={len(parts)}  -> {MANIFEST_OUT}')


def push_to_gcs():
    log(f'Pushing to {GCS_DEST}...')
    cmd = ['gsutil', '-m', 'cp', '-r', str(BASE_OUT), str(MANIFEST_OUT), GCS_DEST + '/']
    log('  ' + ' '.join(cmd))
    subprocess.run(cmd, check=True)
    log('  done.')


def main():
    parts = sorted(BASE_OUT.glob('*.parquet')) if BASE_OUT.exists() else []
    rebuilt = False
    if parts:
        log(f'Base parquet already built at {BASE_OUT} ({len(parts)} part(s)); '
            f'skipping Hail + Polars stages.')
    else:
        build_base_table_hail()
        join_external_scores()
        rebuilt = True

    # ALWAYS regenerate manifest if the parquet was just (re)built. Otherwise
    # only write if missing. This prevents a stale manifest from a prior build
    # being shipped alongside a newer parquet.
    if rebuilt or not MANIFEST_OUT.exists():
        write_manifest()
    else:
        log(f'Manifest already at {MANIFEST_OUT}; not regenerating')

    if '--push' in sys.argv:
        push_to_gcs()
    else:
        log('Skipping GCS push (pass --push to upload).')

    log(f'Done. Inputs staged at {STAGE_DIR}')


if __name__ == '__main__':
    main()
