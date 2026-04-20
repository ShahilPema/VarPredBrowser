#!/usr/bin/env python3
"""Workbench-side: pack curves + figures + k_summary into a single tarball
and push to GCS for retrieval at BCM.

Per-variant fitted_oe parquets are NOT shipped from the workbench. They're
regenerated locally at BCM by `09_apply_curves_local.py` after pulling the
curve pickles. Keeps the AoU egress payload small and free of per-variant
data.
"""
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (
    OUTPUT_DIR, CURVES_DIR, FIG_DIR, GCS_OUTPUTS_PREFIX,
)


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


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
    pack_dir = OUTPUT_DIR / 'pack'
    pack_dir.mkdir(parents=True, exist_ok=True)
    tarball = pack_dir / 'curves_and_figures.tar.zst'

    files = []
    for d in (CURVES_DIR, FIG_DIR):
        for p in sorted(d.rglob('*')):
            if p.is_file():
                files.append(p)
    log(f'Packing {len(files)} files...')

    manifest = {
        'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'files': [
            {'name': str(p.relative_to(OUTPUT_DIR)),
             'size_bytes': p.stat().st_size,
             'md5': md5_file(p)}
            for p in files
        ],
    }
    manifest_path = pack_dir / 'MANIFEST.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    # tar then zstd compress externally (zstd is faster than gzip at similar ratios)
    raw_tar = pack_dir / 'curves_and_figures.tar'
    with tarfile.open(raw_tar, 'w') as tf:
        for p in files:
            tf.add(p, arcname=str(p.relative_to(OUTPUT_DIR)))
        tf.add(manifest_path, arcname='MANIFEST.json')
    log(f'Compressing {raw_tar}...')
    subprocess.run(['zstd', '-3', '--rm', str(raw_tar), '-o', str(tarball)],
                   check=True)

    size_gb = tarball.stat().st_size / 1e9
    log(f'  tarball size: {size_gb:.2f} GB')

    log(f'Pushing {tarball} to {GCS_OUTPUTS_PREFIX}/')
    subprocess.run(['gsutil', '-m', 'cp', str(tarball), f'{GCS_OUTPUTS_PREFIX}/'],
                   check=True)
    log('Done.')


if __name__ == '__main__':
    main()
