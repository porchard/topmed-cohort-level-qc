#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
import numpy as np
import itertools
import seaborn as sns
import textwrap

parser = argparse.ArgumentParser(description='Pivot a table.', add_help=True)
parser.add_argument('--pcs', help='File of PC scores. One column represents a sample ID. Other columns are named PC1, PC2, ...')
parser.add_argument('--variance-explained', dest='variance_explained', help='File of variance explained by each PC.')
parser.add_argument('--covariates', help='QC metrics or other values to be plotted against PC scores. One column represents sample ID. Other column names represent metric names.')
parser.add_argument('--index', help='Name of index column. Must be the same in both files.')
parser.add_argument('--number-pcs', dest='number_pcs', type=int, default=10, help='PCs to plot (default: 10, meaning PCs 1-10 are plotted).')
parser.add_argument('--exclude-covariates', dest='exclude_covariates', nargs='*', default=None, help='List of covariates to exclude from analysis.')
parser.add_argument('--prefix', default='pcs-vs-metrics.', help='Prefix for output files (default: pcs-vs-metrics.)')
args = parser.parse_args()


def convert_str_to_categorical(df):
    """Given a dataframe, convert 'object' columns to pd.Categorical type"""
    categorical_types = df.select_dtypes(include=['object'])
    if len(categorical_types.columns) == 0:
        return df
    else:
        noncategorical_types = df.select_dtypes(exclude=['object'])
        for i in categorical_types.columns:
            categorical_types[i] = pd.Categorical(categorical_types[i], categories=categorical_types[i].unique())
        all_elements = pd.concat([categorical_types, noncategorical_types], axis=1)
        return all_elements.loc[:,df.columns.to_list()]


# given PCs and matrix of covariates, calculate adjusted R^2 and p-values
def variance_in_pcs_explained_by_covariates_handle_missing(pca, covariate_matrix):
    """
    pca is a dataframe with columns = PCs, index = library
    covariate matrix is a dataframe with index = library, columns = covariates
    """
    # verify the indices are the same
    pca = pca.copy().sort_index()
    covariate_matrix = covariate_matrix.copy().sort_index()
    assert(all(covariate_matrix.index == pca.index))
    # run the regressions
    #covariate_matrix = convert_str_to_categorical(covariate_matrix)
    # drop items without any variance
    DROP_COLS = covariate_matrix.nunique().where(lambda x: x<2).dropna().index.to_list()
    covariate_matrix = covariate_matrix.loc[:,~covariate_matrix.columns.to_series().isin(DROP_COLS)]
    results = pd.DataFrame(np.zeros((len(pca.columns), len(covariate_matrix.columns))))
    pvalues = pd.DataFrame(np.ones((len(pca.columns), len(covariate_matrix.columns))), index=pca.columns, columns=covariate_matrix.columns)
    results.index = pca.columns
    results.columns = covariate_matrix.columns
    for PC in pca.columns:
        for var in covariate_matrix.columns:
            predictors = covariate_matrix[var]
            response = pca[PC]
            mask = ~predictors.isnull()
            predictors = predictors[mask]
            response = response[mask]
            if predictors.dtype.name == 'object':
                if predictors.nunique() >= (0.8 * len(predictors)):
                    # skip if there are a large number of categories relative to the number of samples -- as would be the case in sample IDs, for example
                    continue
                predictors = pd.Categorical(predictors, categories=predictors.unique())
                predictors = pd.get_dummies(predictors, drop_first=True)
                predictors.index = response.index
            else:
                predictors = predictors.to_frame()
            # run the OLS
            predictors = sm.add_constant(predictors)
            mod = sm.OLS(response, predictors).fit()
            results.at[PC, var] = mod.rsquared_adj
            pvalues.at[PC, var] = mod.f_pvalue
    return (results, pvalues)



def plot_variance_in_pcs_explained_by_covariates(r2_df, p_df, pcs_variance_explained):
    """
    r2_df is output of variance_in_pcs_explained_by_covariates
    pcs_variance_explained is np array of variance expained by each PC
    """
    r2_df = r2_df.applymap(lambda x: max(x, 0))
    width = 2 + len(r2_df.columns) * 0.3
    height = 2 + len(r2_df) * 0.3

    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(r2_df, vmin=0, vmax=1, cmap='Reds')
    #for r_index in range(len(p_df)):
    #    for c_index in range(len(p_df.columns)):
    #        if p_df.iloc[r_index,c_index] <= 0.05:
    #            ax.text(c_index, r_index, s='*', ha='center', va='center')
    ax.set_xticks(list(range(0, r2_df.shape[1])))
    ax.set_xticklabels(r2_df.columns, fontdict={'rotation': 90})
    ax.set_yticks(list(range(0, len(r2_df))))
    ax.set_yticklabels(['{} ({}%)'.format(i, round(100*j, 2)) for i, j in zip(r2_df.index, pcs_variance_explained)])
    #fraction = 1 if len(r2_df) <= 5 else (5 / len(r2_df))
    fraction = 0.15
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap='Reds'), aspect=10, fraction=fraction)#, pad=0.01)
    #cbar.set_label("max(adj. R^2, 0)")
    cbar.set_label('adj. R^2')
    ax.set_ylim(top=-0.5, bottom=len(r2_df)-0.5)

    return (fig, ax)


pc_scores = pd.read_csv(args.pcs, sep='\t', index_col=0)
pc_var_explained = pd.read_csv(args.variance_explained, sep='\t', index_col=False).set_index('PC').fraction_variance_explained

PCS_PRESENT = list(sorted([[int(i.replace('PC', '')), i] for i in pc_scores.columns], key=lambda x: x[0]))
PCS_KEEP = [i for i in PCS_PRESENT if i[0] <= args.number_pcs]
pc_scores = pc_scores[[i[1] for i in PCS_KEEP]]

covariates = pd.read_csv(args.covariates, sep='\t', index_col=0)
if args.exclude_covariates is not None:
    covariates = covariates.loc[:,~covariates.columns.to_series().isin(args.exclude_covariates)]

IN_BOTH_INDICES = [i for i in pc_scores.index if i in covariates.index]
tmp, p = variance_in_pcs_explained_by_covariates_handle_missing(pc_scores.loc[IN_BOTH_INDICES,:], covariates.loc[IN_BOTH_INDICES,:])
fig, ax = plot_variance_in_pcs_explained_by_covariates(tmp, p, pc_var_explained)
fig.tight_layout()
fig.savefig('{}heatmap.png'.format(args.prefix), bbox_inches='tight')
fig.clf()