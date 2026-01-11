#!/usr/bin/env python3
"""
Test script to validate browser_data.ipynb merge logic.
Run with: python scripts/test_browser_data_merge.py
"""

import hail as hl

# Initialize Hail with minimal resources for testing
hl.init(master='local[4]', quiet=True)

my_bucket = '/storage/zoghbi/home/u235147/merged_vars'

print("=" * 60)
print("STEP 1: Load and inspect base table")
print("=" * 60)

base_ht = hl.read_table(f'{my_bucket}/tmp/constraint_metrics_by_locus_rgc_glaf.ht')
print(f"Base table key: {list(base_ht.key)}")
print(f"Base table row fields: {list(base_ht.row)[:20]}...")  # First 20

# Check what the key fields are
key_fields = list(base_ht.key)
print(f"\nKey fields that CANNOT be in select(): {key_fields}")

print("\n" + "=" * 60)
print("STEP 2: Test column selection")
print("=" * 60)

# Core columns - EXCLUDING key fields
core_cols = ['HGNC', 'chrom', 'pos', 'aa_pos']

# Check if transcript_id exists and is a key or row field
if 'transcript_id' in base_ht.row:
    if 'transcript_id' not in key_fields:
        core_cols.append('transcript_id')
        print("transcript_id is a row field, adding to core_cols")
    else:
        print("transcript_id is a KEY field, NOT adding to core_cols")

# Check region
if 'region' in key_fields:
    print("region is a KEY field, NOT adding to core_cols")
else:
    core_cols.append('region')
    print("region is a row field, adding to core_cols")

# RGC oe/vir columns (excluding key fields)
rgc_oe_vir_cols = [col for col in base_ht.row if ('oe' in col or 'vir' in col) and col not in key_fields]
print(f"Found {len(rgc_oe_vir_cols)} RGC oe/vir columns")

# RGC count columns (excluding key fields)
rgc_count_cols = [col for col in base_ht.row
                  if col.startswith('rgc_')
                  and ('_count' in col or '_obs_' in col or '_prob_mu' in col or '_max_af' in col)
                  and col not in key_fields]
print(f"Found {len(rgc_count_cols)} RGC count/obs columns")

cols_to_keep = core_cols + rgc_oe_vir_cols + rgc_count_cols
print(f"\nTotal columns to select: {len(cols_to_keep)}")
print(f"Columns: {cols_to_keep[:10]}...")

# Verify none are key fields
for col in cols_to_keep:
    if col in key_fields:
        print(f"ERROR: {col} is a key field and should not be in select()!")

print("\n" + "=" * 60)
print("STEP 3: Test select operation")
print("=" * 60)

try:
    test_ht = base_ht.select(*cols_to_keep)
    print("SUCCESS: select() worked!")
    print(f"Result key: {list(test_ht.key)}")
    print(f"Result row fields: {list(test_ht.row)[:10]}...")
except Exception as e:
    print(f"ERROR: {e}")

print("\n" + "=" * 60)
print("STEP 4: Test key_by operation")
print("=" * 60)

try:
    # After select, the key should still be the original key
    # Now re-key to (locus, transcript_id)
    test_ht = test_ht.key_by('locus', 'transcript_id')
    print("SUCCESS: key_by() worked!")
    print(f"New key: {list(test_ht.key)}")
except Exception as e:
    print(f"ERROR: {e}")

print("\n" + "=" * 60)
print("STEP 5: Load and inspect AF tables")
print("=" * 60)

# Load AF tables
rgc_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/no_anno/rgc.ht')
print(f"RGC table key: {list(rgc_ht.key)}")
print(f"RGC has 'alleles': {'alleles' in rgc_ht.row or 'alleles' in rgc_ht.key}")

gnomad_exomes_ht = hl.read_table('/storage/zoghbi/data/sharing/hail_tables/gnomadV4_exomes/gnomadV4exomes_snvs_sex.ht')
print(f"gnomAD exomes key: {list(gnomad_exomes_ht.key)}")

print("\n" + "=" * 60)
print("STEP 6: Test AF aggregation")
print("=" * 60)

try:
    # Test grouping RGC by locus
    rgc_by_locus = rgc_ht.group_by('locus').aggregate(
        rgc_variants = hl.agg.collect(hl.struct(
            alt = rgc_ht.alleles[1],
            af = rgc_ht.info.ALL_AF,
        ))
    )
    print("SUCCESS: RGC group_by/aggregate worked!")
    print(f"Result key: {list(rgc_by_locus.key)}")
    # Show one row
    print("Sample row:")
    rgc_by_locus.show(1)
except Exception as e:
    print(f"ERROR: {e}")

print("\n" + "=" * 60)
print("DONE - All tests completed")
print("=" * 60)
