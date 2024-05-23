"""
perform reconstruction on blind collapsed reads
input blind collapsed reads
output reconstrution result
always recon for anchors
"""

import os
import umap
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from cuml.manifold.umap import UMAP as cuUMAP
from helpers import *

# generate spase matrix from matching, with selection on anchor or target
def get_matrix(match_df, min_a_cnt, max_a_cnt, min_t_cnt, max_t_cnt, anchor, target):
    a_all = match_df.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt')  
    a_sel = a_all.loc[(a_all['total_cnt']>min_a_cnt) & (a_all['total_cnt']<max_a_cnt),]
    t_all = match_df.groupby(target)['cnt'].sum().reset_index(name='total_cnt')  
    t_sel = t_all.loc[(t_all['total_cnt']>min_t_cnt) & (t_all['total_cnt']<max_t_cnt),]
    match_df = match_df[(match_df[anchor].isin(a_sel[anchor])) & (match_df[target].isin(t_sel[target]))]
    a_list = match_df.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt') 
    t_list = match_df.groupby(target)['cnt'].sum().reset_index(name='total_cnt') 
    print('a: {}'.format(len(a_list)))
    print('t: {}'.format(len(t_list)))
    a_dict = dict()
    t_dict = dict()
    for i in range(len(a_list)):
        a_dict[a_list.iloc[i,0]] = i
    for j in range(len(t_list)):
        t_dict[t_list.iloc[j,0]] = j
    a_coo = []
    t_coo = []
    [a_coo.append(a_dict[a]) for a in match_df[anchor]]
    [t_coo.append(t_dict[t]) for t in match_df[target]]
    counts_coo = sp.coo_matrix((match_df['cnt'], (a_coo, t_coo)))
    counts = counts_coo.tocsr().toarray()
    return counts, a_list, t_list

def get_args():
    parser = argparse.ArgumentParser(description='Process recon seq data.')
    parser.add_argument("-f", "--csvpath",
        help="path to the blind_raw_reads_filtered.csv.gz",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-e", "--exptype",
        help="define experiment type (seq or tags).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n", "--n_epochs",
        help="n_epochs paramter",
        type=int,
        default=50000,
    )
    parser.add_argument(
        "-d", "--min_dist",
        help="min_dist paramter",
        type=float,
        default=0.99,
    )
    parser.add_argument(
        "-m", "--min",
        help="I don't know what this does",
        type=int,
        default=0,
    )
    args = parser.parse_args()
    return args

args = get_args()
in_path = args.csvpath
exptype = args.exptype
n_epochs = args.n_epochs
min_dist = args.min_dist
m = args.min
anchor = "V15T" # just for naming
target = "V10T" # just for naming

out_path = f"RUN-n_epochs{n_epochs}-min_dist{min_dist}-min{m}"
if not os.path.exists(out_path):
    os.makedirs(out_path)

print("loading data")
blind_raw = pd.read_csv(os.path.join(in_path,'blind_raw_reads_filtered.csv.gz'))
blind_sum = blind_raw.groupby(['R1_bc', 'R2_bc']).size().reset_index(name='cnt')
if exptype == 'seq':
    blind_sum.columns = [anchor, target, 'cnt'] #if seq
elif exptype == 'tags':
    blind_sum.columns = [target, anchor, 'cnt'] #if tags
del blind_raw

print("plot blind cnt distribution")
a_all = blind_sum.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt')
t_all = blind_sum.groupby(target)['cnt'].sum().reset_index(name='total_cnt')
plot_blind_cnt_distribution(a_all, anchor, path=os.path.join(out_path,anchor+'_blind_cnt_distribution.png'))
plot_blind_cnt_distribution(t_all, target, path=os.path.join(out_path,target+'_blind_cnt_distribution.png'))
del a_all, t_all

print("plot bc covered")
a_cover_bc = blind_sum.groupby(anchor).count()
t_cover_bc = blind_sum.groupby(target).count()
plot_blind_cover_bc_distribution(a_cover_bc, anchor, path=os.path.join(out_path,anchor+'_blind_cover_bc_distribution.png'))
plot_blind_cover_bc_distribution(t_cover_bc, target, path=os.path.join(out_path,target+'_blind_cover_bc_distribution.png'))
del a_cover_bc, t_cover_bc

# generate matrix and reconstruction with CPU
a_min = m
a_max = 1000000
t_min = m
t_max = 1000000
counts, a_sel, t_sel = get_matrix(blind_sum, min_a_cnt=a_min, max_a_cnt=a_max, min_t_cnt=t_min, max_t_cnt=t_max, anchor=anchor, target=target)

reducer = cuUMAP(metric='cosine',
                 n_neighbors=25, 
                 min_dist=min_dist,
                 low_memory=True, 
                 n_components=2, 
                 # random_state=0, 
                 verbose=True, 
                 n_epochs=n_epochs,
                 # output_dens = True,
                 # local_connectivity = 30,
                 learning_rate = 1)
embedding = reducer.fit_transform(np.log1p(counts))

print("output reconstruction result")
a_recon = pd.DataFrame(embedding)
a_recon.columns = ['xcoord','ycoord']
a_recon.insert(loc=0, column=anchor, value=a_sel[anchor])
a_recon.to_csv(os.path.join(out_path,'_recon_loc.csv'), index=False)

print("creating plots")
plot_umap(embedding, path=os.path.join(out_path,'_UMAP.png'))
plot_density(embedding, counts, path=os.path.join(out_path,'_UMAP_density.png'))
plot_convex(embedding, anchor, path=os.path.join(out_path,'_UMAP_convex.png'))
