#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
import seaborn as sns
import re
import logging
import argparse
from adjustText import adjust_text

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Plot PCA, optionally coloring by a chosen variable', add_help=True)
parser.add_argument('--color-by', dest='color_by', type=str, default=None, help='Color by this variable (must be a column in the --metadata file)')
parser.add_argument('--label', type=str, default=None, help='Label libraries listed in this file. TSV, columns: library, label')
parser.add_argument('--id', type=str, default='library', help='This column in all the files represents library')
parser.add_argument('--metadata', default=None, help='File of library metadata. One column must be "library". Required if using --color-by. Includes colnames.')
parser.add_argument('--panel-width', dest='panel_width', type=int, default=5, help='Width of each individual panel (default: 5).')
parser.add_argument('--panel-height', dest='panel_height', type=int, default=5, help='Width of each individual panel (default: 5).')
parser.add_argument('--plot-top', default=12, dest='plot_top', type=int, help='Plot this many PCs. Must be an even number. (default: 12).')
parser.add_argument('pc_scores', help='File of per-library PC scores (columns: library, PC1, PC2, ...). Include colnames.')
parser.add_argument('pc_variance_explained', help='File of variance explained by each PC scores (columns: PC, fraction_variance_explained). Include colnames.')
parser.add_argument('out', help='File to save plot to')
args = parser.parse_args()

if not args.plot_top % 2 == 0:
    raise ValueError('--plot-top must be an even number')

def plot_pca(pc_scores, explained_variance, labels={}, kwargs={}):
    """
    kwargs will be passed through to sns.scatterplot
    """
    number_plots = args.plot_top / 2
    nrows = int(round(np.sqrt(number_plots)))
    ncols = int(np.ceil(number_plots/nrows))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(args.panel_width*ncols, args.panel_height*nrows))
    if number_plots == 1:
        axs = np.array([axs])
    all_texts = []
    for count, ax in enumerate(axs.flatten(), 1):
        first_pc = 'PC{}'.format(2*(count-1)+1)
        second_pc = 'PC{}'.format(2*(count-1)+2)
        xlab = '{} ({}%)'.format(first_pc, round(explained_variance[first_pc]*100, 2))
        ylab = '{} ({}%)'.format(second_pc, round(explained_variance[second_pc]*100, 2))
        ax = sns.scatterplot(x=first_pc, y=second_pc, data=pc_scores, ax=ax, **kwargs)
        if len(labels) > 0:
            labelme = pc_scores[pc_scores.index.isin(labels.keys())]
            texts = [ax.text(x=r[first_pc], y=r[second_pc], s=labels[i]) for i, r in labelme.iterrows()]
            all_texts.append(texts)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    for texts, ax in zip(all_texts, axs.flatten()):
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='red'))
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    for ax in axs.flatten():
        ax.legend().remove()
    fig.tight_layout()
    return (fig, axs)



pc_scores = pd.read_csv(args.pc_scores, sep='\t', low_memory=False, index_col=0)
pc_scores.index = pc_scores.index.to_series().astype(str)
explained_variance = pd.read_csv(args.pc_variance_explained, sep='\t', low_memory=False, index_col=False).set_index('PC').fraction_variance_explained

metadata = None
colors = None
if args.color_by is not None:
    if args.metadata is None:
        raise ValueError('--metadata must be passed if using --color-by')
    metadata = pd.read_csv(args.metadata, sep='\t', index_col=False)
    metadata[args.id] = metadata[args.id].astype(str)
    color_by = metadata.set_index(args.id)[args.color_by]
    colors = color_by[pc_scores.index.to_list()]

labels = dict()
if args.label is not None:
    labels = pd.read_csv(args.label, sep='\t', index_col=False)
    labels[args.id] = labels[args.id].astype(str)
    labels = dict(zip(labels[args.id], labels.label))

fig, axs = plot_pca(pc_scores=pc_scores, explained_variance=explained_variance, labels=labels, kwargs={'alpha': 0.3, 'hue': colors}) if colors is not None else plot_pca(pc_scores=pc_scores, explained_variance=explained_variance, labels=labels, kwargs={'alpha': 0.3})
fig.savefig(args.out, bbox_inches='tight')
fig.clf()
