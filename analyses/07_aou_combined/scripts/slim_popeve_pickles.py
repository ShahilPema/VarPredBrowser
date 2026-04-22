#!/usr/bin/env python3
"""Strip the per-tx `sorted_train_scores` arrays out of curve pickles and
gzip the result.

Both popeve_gp and IRLS-QP (cloglog, quadprog, poisson, bernoulli) fitters
store the per-gene full sorted score universe inside each model dict. The
universe is dataset-independent (same per-gene across all 3 datasets for
a given tag) and trivially rebuildable at apply time from the base score
table. Removing it shrinks each pickle by 10-100x; gzip then trims another
30-50% off the remaining model parameters.

Output filenames end in `.pkl.gz`. The apply paths in 09_apply_curves_local.py
and _popeve_apply.py auto-detect both `.pkl` and `.pkl.gz` (gzip-aware
loader) and detect slim pickles via the missing sorted_train_scores field
to rebuild the universe from scores_df on the fly. Old-style fat pickles
still work unchanged.

By default this script slims every curves_*.pkl in CURVES_DIR. Restrict
with --pattern (e.g. --pattern 'curves_popeve_gp_*.pkl'). --no-gzip writes
plain .pkl instead.

Run on the workbench BEFORE shipping:
    python scripts/slim_popeve_pickles.py
"""
import argparse
import gzip
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import CURVES_DIR  # noqa


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


def slim_one(src: Path, dst_dir: Path, compress: bool = True,
             compresslevel: int = 9) -> Path:
    fits = pickle.load(open(src, 'rb'))
    slim = {tx: {k: v for k, v in d.items() if k != 'sorted_train_scores'}
            for tx, d in fits.items()}
    if compress:
        out = dst_dir / (src.name + '.gz')
        with gzip.open(out, 'wb', compresslevel=compresslevel) as f:
            pickle.dump(slim, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        out = dst_dir / src.name
        with open(out, 'wb') as f:
            pickle.dump(slim, f, protocol=pickle.HIGHEST_PROTOCOL)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src-dir', default=str(CURVES_DIR),
                    help='Directory containing curves_*.pkl')
    ap.add_argument('--dst-dir', default=str(CURVES_DIR / 'slim'),
                    help='Output directory for slim pickles')
    ap.add_argument('--pattern', default='curves_*.pkl',
                    help='Glob pattern (default slims every curves_*.pkl)')
    ap.add_argument('--no-gzip', action='store_true',
                    help='Write plain .pkl instead of .pkl.gz')
    ap.add_argument('--compresslevel', type=int, default=9,
                    help='gzip level 1-9 (default 9, max compression)')
    ap.add_argument('--force', action='store_true',
                    help='Re-slim files that already exist in dst-dir')
    args = ap.parse_args()

    src_dir = Path(args.src_dir)
    dst_dir = Path(args.dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    pkls = sorted(src_dir.glob(args.pattern))
    if not pkls:
        log(f'No pickles matching {args.pattern} in {src_dir}')
        return

    compress = not args.no_gzip
    log(f'Slimming {len(pkls)} pickles -> {dst_dir} '
        f'({"gzip level " + str(args.compresslevel) if compress else "no compression"})')
    total_in = total_out = 0
    skipped = 0
    for src in pkls:
        expected = dst_dir / (src.name + '.gz' if compress else src.name)
        if expected.exists() and not args.force:
            out_mb = expected.stat().st_size / 1024 / 1024
            log(f'  SKIP {src.name} -> {expected.name} exists ({out_mb:.1f} MB)')
            skipped += 1
            continue
        in_mb = src.stat().st_size / 1024 / 1024
        out = slim_one(src, dst_dir, compress=compress,
                       compresslevel=args.compresslevel)
        out_mb = out.stat().st_size / 1024 / 1024
        total_in += in_mb
        total_out += out_mb
        log(f'  {src.name} -> {out.name}: {in_mb:.1f} MB -> {out_mb:.1f} MB '
            f'({100 * out_mb / in_mb:.1f}%)')
    if skipped:
        log(f'Skipped {skipped} file(s) already present in {dst_dir} '
            f'(use --force to re-slim)')

    log(f'Total: {total_in:.0f} MB -> {total_out:.0f} MB '
        f'({100 * total_out / total_in:.1f}%)')


if __name__ == '__main__':
    main()
