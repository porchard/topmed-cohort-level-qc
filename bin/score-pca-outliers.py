#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--normalized-absolute-deviation', dest='normalized_ad', help='')
parser.add_argument('--mahalanobis-distance', dest='md', help='')
parser.add_argument('--mahalanobis-pcs', dest='md_pcs', type=int, help='Number of PCs that were used in calculating mahalanobis distance.')
parser.add_argument('--max-normalized-absolute-deviation', default=5, dest='max_normalized_ad', type=float, help='')
parser.add_argument('--mahalanobis-pvalue-threshold', default=0.01, dest='md_pvalue_threshold', type=float, help='')
parser.add_argument('--prefix', default='qc.', help='Prefix for output files')
args = parser.parse_args()

AD = args.normalized_ad
MD = args.md
PREFIX = args.prefix
MAHALANOBIS_THRESHOLD = args.md_pvalue_threshold
MAX_NORMALIZED_AD = args.max_normalized_ad


md = pd.read_csv(MD, sep='\t', index_col=0)
ad = pd.read_csv(AD, sep='\t', index_col=0)


# plot MD
# https://towardsdatascience.com/multivariate-outlier-detection-in-python-e946cfc843b3
cutoff = chi2.ppf(1-MAHALANOBIS_THRESHOLD, args.md_pcs) ** 0.5
fig, ax = plt.subplots()
sns.histplot(x='mahalanobis_distance', ax=ax, data=md)
ax.set_xlabel('Mahalabonis distance')
ax.set_ylabel('Samples')
ax.axvline(cutoff, color='red', ls='dashed')
fig.tight_layout()
fig.savefig(f'{PREFIX}md-hist.png', dpi=300)
fig.clf()



# plot AD
fig, axs = plt.subplots(ncols=len(ad.columns), figsize=(len(ad.columns)*3, 3))

for col, ax in zip(ad.columns, axs.flatten()):
    sns.histplot(x=col, ax=ax, data=ad)
    ax.set_xlabel('AD / MAD')
    ax.set_ylabel('Samples')
    ax.axvline(MAX_NORMALIZED_AD, color='red', ls='dashed')
fig.tight_layout()
fig.savefig(f'{PREFIX}ad-hist.png', dpi=300)
fig.clf()



assert(all(ad.index == md.index))
all_samples = set(ad.index)
fail_md = set(md[md.pvalue<MAHALANOBIS_THRESHOLD].index)
fail_ad = set(ad[ad.max(axis=1)>MAX_NORMALIZED_AD].index)
pca_outliers = fail_md.union(fail_ad)
nonoutliers = all_samples.difference(pca_outliers)

outlier_status = [[i, 'outlier'] for i in pca_outliers] + [[i, 'nonoutlier'] for i in nonoutliers]
outlier_status = pd.DataFrame(outlier_status).to_csv(f'{PREFIX}outlier-status.txt', sep='\t', index=False, header=False)