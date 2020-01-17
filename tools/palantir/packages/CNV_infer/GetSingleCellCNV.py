# load modules
import sys
import pickle
import json
import os
import numpy as np
import pandas as pd
import pickle
from collections import defaultdict, Counter, OrderedDict
from imp import reload
from itertools import combinations
from copy import deepcopy
import matplotlib.pyplot as plt

from plotCnvArr import plotcnv
from inferCNVFast import infer_cnv

# Plotting imports
import matplotlib
import seaborn as sns
import matplotlib.gridspec as gridspec

import re
import seqc
import harmony
import palantir
import phenograph

out_dir=sys.argv[1]
counts_file=sys.argv[2]
counts_g_file=sys.argv[3]
counts_bc_file=sys.argv[4]
cnv_file=sys.argv[5]
plot_file=sys.argv[6]

sns.set_style('white')

#cnv_df3 = pd.read_hdf(out_dir + 'cnv_df.epithelial.wo_dying_cells.070519.h5',key='cnv_df')

from scipy.io import mmread
fn = out_dir + counts_file
counts = mmread(fn)

fn = out_dir + counts_g_file
with open(fn, 'r') as f:
    g = [i.strip() for i in f.readlines()]

fn = out_dir + counts_bc_file
with open(fn, 'r') as f:
    bc = [i.strip() for i in f.readlines()]

counts = pd.DataFrame(counts.tocsr().toarray(), index=bc, columns = g)

#normal_df = pd.read_hdf('/data/peer/chanj3/HTA.combined.052619/SEQC/LUNG_ADENOCARCINOMA_CLEAN.h5',key='DF_EPITHELIAL_NOR_TUMOR')
#normal_df = normal_df.loc[normal_df.index.get_level_values('Legend').str.contains('NORMAL')]
#meta_normal = normal_df.index
#normal_df.index = [i + '_' + j for i,j in zip(meta_normal.get_level_values('Legend'), meta_normal.get_level_values('Cell ID'))]
#counts = pd.concat([counts,normal_df],axis=0).fillna(0)

import scanpy as sc
adata = sc.AnnData(X = counts)
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_per_cell(adata) #counts_per_cell_after= np.median(adata.obs['original_total_counts']))
sc.pp.filter_genes(adata, min_cells=10)
norm_unfiltered_df = pd.DataFrame(adata.X, index=adata.obs_names, columns = adata.var_names)

data_df = np.log2(norm_unfiltered_df+0.1)

import pickle

batch = [re.sub('_[0-9]+$','',i) for i in data_df.index]
batch = pd.Series(batch, index = data_df.index)

genePosFile = '/home/chanj3/data/HTA.NSCLC_epithelial.plasticity.070219/notebooks/infercnv/genePos.wo_MT.txt'
chrFile = '/home/chanj3/data/HTA.NSCLC_epithelial.plasticity.070219/notebooks/infercnv/chrNameLength.txt'

cnv_df = infer_cnv(genePosFile,chrFile,data_df, narounds=100)

#####MUST CHANGE FOR GETTING NORMAL/WILDTYPE
diploid_mean = cnv_df.loc[cnv_df.index.str.contains('_N_|_LN_minus'),:].mean(axis=0)
diploid_std = cnv_df.loc[cnv_df.index.str.contains('_N_|_LN_minus'),:].std(axis=0)

cnv_df2 = cnv_df.subtract(diploid_mean, axis=1).div(diploid_std, axis=1)
    

allv=np.ravel(cnv_df2)
allv=allv[np.nonzero(allv)]
maxv=np.percentile(allv,99)
minv=np.percentile(allv,1)
print(maxv, minv)

cnv_df2 = cnv_df2.sub(cnv_df2.median(axis=1), axis=0)    
wt_std = cnv_df2.loc[cnv_df2.index.str.contains('_N_|_LN_minus'),:].std(axis=0)

from copy import deepcopy
tmp = deepcopy(cnv_df2.values)
tmp[np.where(np.abs(cnv_df2).lt(1.5*wt_std, axis=1))] = 0
cnv_df3 = pd.DataFrame(tmp, index=cnv_df2.index, columns=cnv_df2.columns)

cnv_df3.shape

store = pd.HDFStore(out_dir + cnv_file)
store['cnv_df'] = cnv_df3.astype(np.float32)
store.close()

pal_tn = sns.color_palette('Set1', 2)
lut_tn = dict(zip([False, True], pal_tn))

pal_tp = sns.color_palette('Set1', len(set(batch)))
lut_tp = dict(zip(sorted(list(set(batch)), key = lambda x: '_N_' in x or 'minus' in x), pal_tp))

pal_chr = sns.color_palette('deep', 23)
lut_chr = dict(zip(list(range(1,24)), pal_chr))

row_colors1 = [lut_tp[ re.sub('_[0-9]+$','',i) ] for i in cnv_df3.index]
row_colors2 = [lut_tn[ 'NORMAL' in i ] for i in cnv_df3.index]
col_colors = [lut_chr[i] for i in cnv_df3.columns.get_level_values(0)]

def normalize(data):
    means = np.ravel(data.mean(axis=1))
    std = np.ravel(data.std(axis=1))
    ret = ((data.T - means) / std).T
    return ret

g = sns.clustermap(cnv_df3, col_cluster = False, center = 0, vmin=-4, vmax=4, method='ward',
                   cmap = plt.cm.seismic,
                   col_colors = col_colors,
                   row_colors = [row_colors1,row_colors2])
ax = g.ax_heatmap
ax.set_xlabel("")
ax.set_xticklabels("")

fname = out_dir + plot_file
g.savefig(fname, dpi=300)
