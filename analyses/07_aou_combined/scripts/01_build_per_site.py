#!/usr/bin/env python3
"""Workbench-side: build per_site_{gnomad_only,aou_only,combined}.parquet.

Reads the consolidated base_table.parquet (already has gnomAD AC + PASS, scores,
gnomAD-derived position filter flags, roulette MR, ENST/HGNC). Joins AoU per-site
AC + PASS + AoU-recomputed coverage/AB flags from the AoU VAT/coverage MT, then
writes three per-site parquets following the canonical hardfilter+cov+ab scenario
shape from `02_poisson_quadprog_4models/scripts/01_extract_per_site.py`.

For each dataset:
  - apply observed rule + position filters
  - fit Poisson scaling factor k_auto and k_chrX on synonymous scaling sites
  - expected = 1 - exp(-k * roulette_AR_MR), per region
  - keep missense rows only
  - write (locus_str, alleles_str, transcript, gene_symbol, observed, expected)

Requires Hail (for the AoU VAT join) — run on a dataproc cluster.
After this step everything is pure-python.
"""
import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (
    BASE_TABLE_PARQUET, PER_SITE_DIR, DATASETS, AOU_CDR_VERSION, per_site_path,
)


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# -----------------------------------------------------------------------------
# AoU VAT join — workbench-specific
# -----------------------------------------------------------------------------
# The AoU Variant Annotation Table and Coverage MatrixTable live at
# workspace-specific paths. Set these via env vars on the workbench:
#   AOU_VAT_HT  — Hail Table with (locus, alleles, ac, filters) for AoU
#   AOU_COV_MT  — Hail MatrixTable with per-sample coverage for AoU
# If the user already has aggregated coverage/AB stats, they can drop in a
# parquet at AOU_COVAB_PARQUET with (locus_str, f_lowcov_aou, f_ab_low_aou).

AOU_VAT_HT = os.environ.get('AOU_VAT_HT', '')
AOU_COV_MT = os.environ.get('AOU_COV_MT', '')
AOU_COVAB_PARQUET = os.environ.get('AOU_COVAB_PARQUET', '')

LOWCOV_DP_THRESHOLD = 20.0   # mirrors gnomAD's f_lowcov definition
AB_LOW_THRESHOLD = 0.3       # mirrors gnomAD's f_ab_low definition

SPARK_TMP = os.environ.get('SPARK_TMP_DIR', '/tmp/aou_combined_spark')


def init_hail():
    import hail as hl
    cpus = int(os.environ.get('HAIL_CPUS', '80'))  # AoU workbench cap is 96 vCPU
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


def load_base_with_aou(hl):
    """Read shipped base parquet via Spark, join AoU AC + PASS + cov/AB flags."""
    log(f'Loading base parquet: {BASE_TABLE_PARQUET}')
    base_df = hl.import_table(
        f'{BASE_TABLE_PARQUET}/*.parquet',
        no_header=False, force=False, types={},
        impute=False,
    )
    # Importing parquet via Hail's spark interop:
    from pyspark.sql import SparkSession
    spark = SparkSession.builder.getOrCreate()
    sdf = spark.read.parquet(str(BASE_TABLE_PARQUET) + '/*.parquet')
    base = hl.Table.from_spark(sdf)
    # locus_str + alleles_str -> Hail locus/alleles
    base = base.annotate(
        _contig=base.locus_str.split(':')[0],
        _pos=hl.int32(base.locus_str.split(':')[1]),
    )
    base = base.annotate(
        locus=hl.locus(base._contig, base._pos, reference_genome='GRCh38'),
        alleles=base.alleles_str.replace('\\[', '').replace('\\]', '')
                                .replace('"', '').split(','),
    )
    base = base.key_by('locus', 'alleles').drop('_contig', '_pos')
    log(f'  base rows (post-key): {base.count():,}')

    if AOU_VAT_HT:
        log(f'Joining AoU VAT: {AOU_VAT_HT}')
        vat = hl.read_table(AOU_VAT_HT)
        if 'locus' not in vat.key:
            vat = vat.key_by('locus', 'alleles')
        _v = vat[base.locus, base.alleles]
        # AoU schema typically has total AC and a filters set/PASS bool — adapt as needed
        base = base.annotate(
            aou_in_vat=hl.is_defined(_v),
            aou_ac=hl.or_else(_v.AC, 0),                    # adjust field name to AoU schema
            aou_is_pass=hl.or_else(hl.len(_v.filters) == 0, False),
        )
    else:
        log('  AOU_VAT_HT not set; aou fields will be empty (gnomad-only leg still works)')
        base = base.annotate(
            aou_in_vat=False, aou_ac=hl.int32(0), aou_is_pass=False,
        )

    if AOU_COVAB_PARQUET and Path(AOU_COVAB_PARQUET).exists():
        log(f'Joining pre-computed AoU cov/AB flags: {AOU_COVAB_PARQUET}')
        sdf2 = spark.read.parquet(AOU_COVAB_PARQUET)
        cov_ht = hl.Table.from_spark(sdf2)
        if 'locus' not in cov_ht.key:
            cov_ht = cov_ht.key_by('locus_str')
        _c = cov_ht[base.locus_str]
        base = base.annotate(
            f_lowcov_aou=hl.or_else(_c.f_lowcov_aou, False),
            f_ab_low_aou=hl.or_else(_c.f_ab_low_aou, False),
        )
    elif AOU_COV_MT:
        log(f'Computing AoU cov/AB flags from {AOU_COV_MT}')
        cov_mt = hl.read_matrix_table(AOU_COV_MT)
        # Aggregate per-locus median DP and median het-call AB
        agg = cov_mt.annotate_rows(
            median_dp=hl.agg.approx_median(cov_mt.DP),
            median_ab=hl.agg.filter(
                cov_mt.GT.is_het(),
                hl.agg.approx_median(cov_mt.AD[1] / hl.sum(cov_mt.AD)),
            ),
        ).rows()
        agg = agg.annotate(
            f_lowcov_aou=hl.or_else(agg.median_dp < LOWCOV_DP_THRESHOLD, False),
            f_ab_low_aou=hl.or_else(agg.median_ab < AB_LOW_THRESHOLD, False),
        ).key_by('locus')
        _c = agg[base.locus]
        base = base.annotate(
            f_lowcov_aou=hl.or_else(_c.f_lowcov_aou, False),
            f_ab_low_aou=hl.or_else(_c.f_ab_low_aou, False),
        )
    else:
        log('  No AoU coverage source set; f_lowcov_aou/f_ab_low_aou = False')
        base = base.annotate(f_lowcov_aou=False, f_ab_low_aou=False)

    return base


# -----------------------------------------------------------------------------
# Poisson k scaling — verbatim from 02_poisson_quadprog_4models extract
# -----------------------------------------------------------------------------

def compute_scaling_factor_hail(hl, ht, obs_field, chrx_region,
                                raw_mr_field='roulette_AR_MR',
                                scaling_flag='roulette_syn_scaling_site'):
    if chrx_region == 'chrX_nonpar':
        scaling_ht = ht.filter(
            ht.is_chrX_nonPAR &
            (ht.most_deleterious_consequence_cds == 'synonymous_variant') &
            hl.is_defined(ht[raw_mr_field]) &
            (ht[raw_mr_field] > 0)
        )
    else:
        scaling_ht = ht.filter(
            (ht[scaling_flag] == True) &
            hl.is_defined(ht[raw_mr_field]) &
            (ht[raw_mr_field] > 0)
        )
        if chrx_region == 'autosomes_par':
            scaling_ht = scaling_ht.filter(~scaling_ht.is_chrX_nonPAR)

    agg_ht = scaling_ht.group_by(raw_mr_field).aggregate(
        polymorphic=hl.agg.sum(scaling_ht[obs_field]),
        total=hl.agg.count(),
    )
    agg_df = agg_ht.to_pandas()

    rates = agg_df[raw_mr_field].values
    totals = agg_df['total'].values
    observed = float(agg_df['polymorphic'].sum())

    def objective(x):
        expected = np.sum((1 - np.exp(-x[0] * rates)) * totals)
        return abs(expected - observed)

    res = minimize(objective, [1.0], method='Nelder-Mead',
                   options={'xatol': 1e-3, 'maxiter': 10000})
    return float(res.x[0])


# -----------------------------------------------------------------------------
# Per-dataset extraction
# -----------------------------------------------------------------------------

def extract_one(hl, dataset, base):
    out_path = per_site_path(dataset)
    if (out_path / '_SUCCESS').exists():
        log(f'SKIP {dataset}: {out_path} exists')
        return

    log(f'=== {dataset} ===')
    ht = base

    if dataset == 'gnomad_only':
        ht = ht.filter(~(ht.gnomad_in_table & (ht.gnomad_ac == 0)))
        ht = ht.annotate(observed=hl.if_else(ht.gnomad_ac > 0, 1, 0))
        ht = ht.filter(~ht.f_lowcov)
        ht = ht.filter(~ht.f_ab_low)

    elif dataset == 'aou_only':
        ht = ht.filter(~(ht.aou_in_vat & (ht.aou_ac == 0)))
        ht = ht.annotate(observed=hl.if_else(ht.aou_ac > 0, 1, 0))
        ht = ht.filter(~ht.f_lowcov_aou)
        ht = ht.filter(~ht.f_ab_low_aou)

    elif dataset == 'combined':
        # Drop site only if BOTH cohorts called it AND BOTH failed hardfilter
        both_called_both_fail = (
            ht.gnomad_in_table & (ht.gnomad_ac > 0) & ~ht.gnomad_is_pass &
            ht.aou_in_vat & (ht.aou_ac > 0) & ~ht.aou_is_pass
        )
        ht = ht.filter(~both_called_both_fail)
        ht = ht.annotate(observed=hl.if_else(
            ((ht.gnomad_ac > 0) & ht.gnomad_is_pass) |
            ((ht.aou_ac > 0) & ht.aou_is_pass),
            1, 0,
        ))
        # Position filter: drop only if both cohorts have bad coverage / bad AB
        ht = ht.filter(~(ht.f_lowcov & ht.f_lowcov_aou))
        ht = ht.filter(~(ht.f_ab_low & ht.f_ab_low_aou))

    else:
        raise ValueError(dataset)

    log('  fitting k_auto / k_chrX on synonymous sites...')
    k_auto = compute_scaling_factor_hail(hl, ht, 'observed', 'autosomes_par')
    k_chrX = compute_scaling_factor_hail(hl, ht, 'observed', 'chrX_nonpar')
    log(f'  k_auto={k_auto:.6f}  k_chrX={k_chrX:.6f}')

    ht = ht.annotate(expected=hl.if_else(
        ht.is_chrX_nonPAR,
        1.0 - hl.exp(-1.0 * ht.roulette_AR_MR * k_chrX),
        1.0 - hl.exp(-1.0 * ht.roulette_AR_MR * k_auto),
    ))

    ht = ht.filter(
        (ht.most_deleterious_consequence_cds == 'missense_variant') &
        hl.is_defined(ht.expected) & (ht.expected > 0)
    )

    out_cols = ['locus_str', 'alleles_str', 'transcript', 'gene_symbol',
                'observed', 'expected']
    ht.key_by().select(*out_cols).to_spark().write.mode('overwrite').parquet(str(out_path))
    log(f'  wrote {out_path}  ({ht.count():,} rows)')


def write_synonymous_slice(hl, base):
    """Write a synonymous-only slice with the columns 02_diagnostic_plots needs."""
    syn_out = PER_SITE_DIR / 'synonymous_augmented.parquet'
    if (syn_out / '_SUCCESS').exists():
        log(f'SKIP synonymous slice ({syn_out} exists)')
        return
    log('Building synonymous_augmented.parquet for diagnostics...')
    syn = base.filter(
        base.most_deleterious_consequence_cds == 'synonymous_variant'
    ).filter(hl.is_defined(base.roulette_AR_MR) & (base.roulette_AR_MR > 0))
    cols = ['locus_str', 'transcript', 'roulette_AR_MR', 'is_chrX_nonPAR',
            'roulette_syn_scaling_site',
            'gnomad_in_table', 'gnomad_ac', 'gnomad_is_pass',
            'aou_in_vat', 'aou_ac', 'aou_is_pass',
            'f_lowcov', 'f_ab_low', 'f_lowcov_aou', 'f_ab_low_aou']
    syn.key_by().select(*cols).to_spark().write.mode('overwrite').parquet(str(syn_out))
    log(f'  wrote {syn_out}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', default=','.join(DATASETS),
                    help=f'comma-separated subset of {DATASETS}')
    args = ap.parse_args()

    requested = args.datasets.split(',')
    for d in requested:
        if d not in DATASETS:
            raise ValueError(f'unknown dataset {d}; expected {DATASETS}')

    log(f'AoU CDR pin: {AOU_CDR_VERSION}')
    os.makedirs(SPARK_TMP, exist_ok=True)

    hl = init_hail()
    base = load_base_with_aou(hl)
    # checkpoint augmented base once
    ckpt = Path(SPARK_TMP) / 'base_with_aou.ht'
    if not (ckpt / '_SUCCESS').exists():
        log(f'Checkpointing augmented base -> {ckpt}')
        base = base.checkpoint(str(ckpt), overwrite=True)
    else:
        base = hl.read_table(str(ckpt))
        log(f'Reusing augmented base from {ckpt}')

    write_synonymous_slice(hl, base)

    for d in requested:
        extract_one(hl, d, base)

    hl.stop()
    log(f'Done. Per-site parquets under {PER_SITE_DIR}')


if __name__ == '__main__':
    main()
