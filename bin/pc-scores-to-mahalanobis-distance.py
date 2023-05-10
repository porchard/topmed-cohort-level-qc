#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import sys
from scipy.spatial import distance
from scipy.stats import chi2
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('pc_scores', help='PC score file')
parser.add_argument('--use-pcs', default=5, dest='use_pcs', type=int, help='Use this many PCs for calculating mahalanobis distance')
args = parser.parse_args()

PC_SCORES = args.pc_scores
USE_PCS = args.use_pcs


def mahalanobis_distance(df):
    """
    Input: df (rows = samples, columns = variables)
    TODO: test
    """
    center = df.mean()
    # get inverse of covariance matrix
    covariance  = np.cov(df, rowvar=False)
    covariance_inverse = np.linalg.inv(covariance)
    dist = df.apply(lambda x: distance.mahalanobis(x, center, covariance_inverse), axis=1)
    return dist


pc_scores = pd.read_csv(PC_SCORES, sep='\t', index_col=0).loc[:,[f'PC{i}' for i in range(1, USE_PCS+1)]]

dist = mahalanobis_distance(pc_scores).rename('mahalanobis_distance').to_frame()
dist['pvalue'] = 1-chi2.cdf(np.power(dist.mahalanobis_distance, 2), USE_PCS)
dist.to_csv(sys.stdout, sep='\t', index=True)
