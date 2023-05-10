#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('pc_scores', help='PC score file')
parser.add_argument('--normalize-by-mad', dest='normalize_by_mad', default=False, action='store_true', help='Divide by MAD')
args = parser.parse_args()

PC_SCORES = args.pc_scores

pc_scores = pd.read_csv(PC_SCORES, sep='\t', index_col=0)

ad = pc_scores.transform(lambda x: (x-x.median()).abs())
if not args.normalize_by_mad:
    ad.to_csv(sys.stdout, sep='\t', index=True)
else:
    mad = ad.median()
    ad_over_mad = ad / mad
    ad_over_mad.to_csv(sys.stdout, sep='\t', index=True)

