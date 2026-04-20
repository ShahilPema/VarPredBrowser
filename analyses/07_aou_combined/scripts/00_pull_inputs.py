#!/usr/bin/env python3
"""Workbench-side: pull base_table.parquet + MANIFEST.json from GCS, verify md5s.

Run on the AoU workbench once after `git clone` and before `01_build_per_site.py`.
"""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import INPUTS_DIR, BASE_TABLE_PARQUET, INPUT_MANIFEST, GCS_INPUTS_PREFIX


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def md5_file(path: Path, chunk=8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    log(f'Pulling from {GCS_INPUTS_PREFIX} -> {INPUTS_DIR}')
    subprocess.run(
        ['gsutil', '-m', 'cp', '-r',
         f'{GCS_INPUTS_PREFIX}/base_table.parquet',
         f'{GCS_INPUTS_PREFIX}/MANIFEST.json',
         str(INPUTS_DIR) + '/'],
        check=True,
    )

    log(f'Verifying via {INPUT_MANIFEST}...')
    with open(INPUT_MANIFEST) as f:
        manifest = json.load(f)

    bad = 0
    for entry in manifest['partitions']:
        local = BASE_TABLE_PARQUET / entry['name']
        if not local.exists():
            log(f'  MISSING: {local}')
            bad += 1
            continue
        if local.stat().st_size != entry['size_bytes']:
            log(f'  SIZE MISMATCH: {local}')
            bad += 1
            continue
        actual = md5_file(local)
        if actual != entry['md5']:
            log(f'  MD5 MISMATCH: {local} got={actual} expected={entry["md5"]}')
            bad += 1

    if bad:
        log(f'{bad} partitions failed verification.')
        sys.exit(1)
    log(f'All {len(manifest["partitions"])} partitions verified '
        f'({manifest["n_rows"]:,} rows).')


if __name__ == '__main__':
    main()
