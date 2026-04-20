#!/usr/bin/env python3
"""Bernoulli cloglog IRLS-QP fits across {gnomad_only, aou_only, combined} x 8 tags.

Outputs: output/curves/curves_cloglog_<tag>_<dataset>.pkl
"""
from _fit_common import standard_main

if __name__ == '__main__':
    standard_main('cloglog')
