#!/usr/bin/env python3
"""Strip the per-tx `sorted_train_scores` arrays out of popeve_gp pickles.

The universe is dataset-independent (same per-gene across all 3 datasets for
a given tag) and trivially rebuildable at apply time from the base score
table. Removing it shrinks each pickle from ~600 MB to ~10-30 MB so they
clear AoU's 100 MB single-file egress limit without splitting.

The apply path (_popeve_apply.py:apply_curves) detects slim pickles via
the missing `sorted_train_scores` field and rebuilds the universe from
scores_df on the fly. Old-style pickles still work unchanged.

Run on the workbench BEFORE tarring/pushing:
    python scripts/slim_popeve_pickles.py
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import CURVES_DIR  # noqa


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


def slim_one(src: Path, dst_dir: Path) -> Path:
    fits = pickle.load(open(src, 'rb'))
    slim = {}
    for tx, d in fits.items():
        new = {k: v for k, v in d.items() if k != 'sorted_train_scores'}
        slim[tx] = new
    out = dst_dir / src.name
    with open(out, 'wb') as f:
        pickle.dump(slim, f, protocol=pickle.HIGHEST_PROTOCOL)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src-dir', default=str(CURVES_DIR),
                    help='Directory containing curves_popeve_gp_*.pkl')
    ap.add_argument('--dst-dir', default=str(CURVES_DIR / 'slim'),
                    help='Output directory for slim pickles')
    args = ap.parse_args()

    src_dir = Path(args.src_dir)
    dst_dir = Path(args.dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    pkls = sorted(src_dir.glob('curves_popeve_gp_*.pkl'))
    if not pkls:
        log(f'No popeve_gp pickles in {src_dir}')
        return

    log(f'Slimming {len(pkls)} pickles -> {dst_dir}')
    total_in = total_out = 0
    for src in pkls:
        in_mb = src.stat().st_size / 1024 / 1024
        out = slim_one(src, dst_dir)
        out_mb = out.stat().st_size / 1024 / 1024
        total_in += in_mb
        total_out += out_mb
        log(f'  {src.name}: {in_mb:.1f} MB -> {out_mb:.1f} MB '
            f'({100 * out_mb / in_mb:.1f}%)')

    log(f'Total: {total_in:.0f} MB -> {total_out:.0f} MB '
        f'({100 * total_out / total_in:.1f}%)')


if __name__ == '__main__':
    main()
