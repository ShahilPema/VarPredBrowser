#!/usr/bin/env python3
"""PLM-alone Bernoulli (no offset) baseline across 3 datasets x 8 tags.

This is `bernoulli_no_offset_irls_qp` — fits P(observed | score) without using
the per-site mutation rate. Useful as a sanity baseline against the offset
methods (cloglog, poisson, quadprog).

Outputs: output/curves/curves_plm_alone_<tag>_<dataset>.pkl
"""
from _fit_common import standard_main

if __name__ == '__main__':
    standard_main('plm_alone')
