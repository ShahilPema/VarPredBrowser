#!/usr/bin/env python3
"""Poisson IRLS-QP fits across 3 datasets x 4 in-house tags.

Outputs: output/curves/curves_poisson_<tag>_<dataset>.pkl
"""
from _fit_common import standard_main

if __name__ == '__main__':
    standard_main('poisson')
