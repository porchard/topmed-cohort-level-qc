#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

OUTPUT, RNASEQC = sys.argv[1:3]
OUTLIER_STATUS = sys.argv[3:]


def plot_on_log_scale(x):
    # determine how many orders of magnitude the data ranges
    if x.min() <= 0:
        return False
    tmp = np.log10(x)
    if (tmp.max() - tmp.min()) >= 2:
        return True
    else:
        return False


rnaseqc = pd.read_csv(RNASEQC, sep='\t', index_col=0)
samples = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['tor', 'outlier_status', 'group']).assign(group=os.path.basename(f).replace('.outlier-status.txt', '')) for f in OUTLIER_STATUS])
rnaseqc = rnaseqc.melt(ignore_index=False).reset_index().merge(samples)

NCOLS = 5
N_METRICS = rnaseqc.variable.nunique()
NROWS = divmod(N_METRICS, NCOLS)[0] + int(divmod(N_METRICS, NCOLS)[1]>0)
fig, axs = plt.subplots(nrows=NROWS, ncols=NCOLS, figsize=(NCOLS*6, NROWS*4))
for metric, ax in itertools.zip_longest(rnaseqc.variable.unique(), axs.flatten()):
    if metric is None:
        ax.remove()
        continue
    df = rnaseqc[rnaseqc.variable==metric]
    
    if len(df) == 0:
        continue
    sns.boxplot(x='group', y='value', data=df, color='grey', showfliers=False, ax=ax)
    sns.stripplot(x='group', y='value', data=df, hue='outlier_status', alpha=0.1, ax=ax, palette={'outlier': 'red', 'nonoutlier': 'black'})
    ax.set_ylim(bottom=0.8*df.value.min(), top=1.2*df.value.max())
    ax.legend().remove()
    if plot_on_log_scale(df.value):
        ax.set_yscale('log')
    for t in ax.get_xticklabels():
        t.set(rotation=45, ha='right')
    ax.set_xlabel('')
    ax.set_ylabel(metric)
    ax.set_title(metric)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=300, facecolor='white')
