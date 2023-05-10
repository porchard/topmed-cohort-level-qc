#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from qtl import norm
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Perform PCA on a gene count matrix', add_help=True)
parser.add_argument('mat', help='Matrix (rows = genes, columns = samples) in parquet format.')
parser.add_argument('--sample-list', default=None, help='File listing samples to include (if not given, use all)')
args = parser.parse_args()

# for any 0s: fill with the (non-zero minimum of that feature) / 2
def fill_zeros(y):
    min_nonzero_val = min(y[y > 0])
    return np.maximum(y, min_nonzero_val/2)


logging.info('Reading {}'.format(args.mat))
samples = pd.read_csv(args.sample_list, header=None)[0].to_list() if args.sample_list is not None else None
mat = pd.read_parquet(args.mat, columns=samples)
logging.info('Finished reading {}'.format(args.mat))

logging.info('Performing DESeq2 normalization')
mat = norm.deseq2_normalized_counts(mat)

logging.info('Dropping genes w/ only 0s')
mat = mat[mat.max(axis=1)>0]

logging.info('Removing lowly expressed genes, filling 0s, and log10-transforming')
mat = mat.T
mask = (mat >= 10).mean() > 0.1
mat = mat.transform(fill_zeros)
mat = np.log10(mat.loc[:,mask])


logging.info('Performing PCA')
logging.info('Centering matrix') # although I do believe the PCA() centers (but does not scale) automatically...
mat = mat - mat.mean()
pca = PCA().fit(mat)
logging.info('Finished PCA.')
explained_variance = pca.explained_variance_ratio_

logging.info('Outputting PC scores to file PC-scores.txt')
transformed = pd.DataFrame(pca.transform(mat), columns=['PC{}'.format(i) for i in range(1, len(explained_variance)+1)])
transformed.index = mat.index
transformed.reset_index().to_csv('PC-scores.txt', sep='\t', index=False)

logging.info('Outputting PC variance explained to file PC-variance-explained.txt')
explained_variance = pd.DataFrame([['PC'+str(count), e] for count, e in enumerate(explained_variance, 1)], columns=['PC', 'fraction_variance_explained'])
explained_variance.to_csv('PC-variance-explained.txt', sep='\t', index=False)

logging.info('Done.')
