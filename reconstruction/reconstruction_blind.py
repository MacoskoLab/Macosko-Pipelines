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
from helpers import *

# generate spase matrix from matching, with selection on anchor or target
def get_matrix(match_df, min_a_cnt, max_a_cnt, min_t_cnt, max_t_cnt, anchor, target):
    a_all = match_df.groupby(anchor)['cnt'].size().reset_index(name='bead_cnt')   
    a_sel = a_all.loc[(a_all['bead_cnt']>min_a_cnt) & (a_all['bead_cnt']<max_a_cnt),]
    t_all = match_df.groupby(target)['cnt'].size().reset_index(name='bead_cnt')  
    t_sel = t_all.loc[(t_all['bead_cnt']>min_t_cnt) & (t_all['bead_cnt']<max_t_cnt),]
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
    counts = counts_coo.tocsr()
    return counts, a_list, t_list

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir",
        help="input data folder",
        type=str,
        required=True,
    )
    parser.add_argument("-o", "--out_dir",
        help="output data folder",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-e", "--exptype",
        help="define experiment type (seq or tags)",
        type=str,
        required=True,
    )
    # Optional args
    parser.add_argument(
        "-n", "--n_neighbors",
        help="the number of neighboring sample points used for manifold approximation",
        type=int,
        default=25,
    )
    parser.add_argument(
        "-m", "--metric",
        help="the metric to use to compute distances in high dimensional space",
        type=str,
        default="cosine",
    )
    parser.add_argument(
        "-N", "--n_epochs",
        help="the number of training epochs to be used in optimizing the low dimensional embedding",
        type=int,
        default=50000,
    )
    parser.add_argument(
        "-d", "--min_dist",
        help="the effective minimum distance between embedded points",
        type=float,
        default=0.99,
    )
    parser.add_argument(
        "-c", "--core",
        help="define core type to use (CPU or GPU)",
        type=str,
        default="CPU",
    )
    args = parser.parse_args()
    return args

args = get_args()
in_dir = args.in_dir
out_dir = args.out_dir
out_dir = os.path.join(out_dir,f'N{args.n_epochs}_n{args.n_neighbors}_d{min_dist}')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
exptype = args.exptype
anchor = "anchor"
target = "target"

print("Loading data...")
blind_raw = pd.read_csv(os.path.join(in_dir, 'blind_raw_reads_filtered.csv.gz'))
blind_sum = blind_raw.groupby(['R1_bc', 'R2_bc']).size().reset_index(name='cnt')
if exptype == 'seq':
    blind_sum.columns = [anchor, target, 'cnt']
elif exptype == 'tags':
    blind_sum.columns = [target, anchor, 'cnt']
else:
    assert False, f"unknown exptype ({exptype})"

print("Plotting blind cnt distribution...")
a_all = blind_sum.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt')
t_all = blind_sum.groupby(target)['cnt'].sum().reset_index(name='total_cnt')
blind_cnt_distribution(a_all, anchor, os.path.join(out_dir, f'blind_cnt_distribution_{anchor}.png'))
blind_cnt_distribution(t_all, target, os.path.join(out_dir, f'blind_cnt_distribution_{target}.png'))

print("Plotting bc covered...")
a_cover_bc = blind_sum.groupby(anchor).count()
t_cover_bc = blind_sum.groupby(target).count()
blind_cover_bc_distribution(a_cover_bc, anchor, os.path.join(out_dir, f'blind_cover_bc_distribution_{anchor}.png'))
blind_cover_bc_distribution(t_cover_bc, target, os.path.join(out_dir, f'blind_cover_bc_distribution_{target}.png'))

print("Generating matrix...")
a_min = 0
a_max = 1000
t_min = 0
t_max = 1000
counts, a_sel, t_sel = get_matrix(blind_sum, min_a_cnt=a_min, max_a_cnt=a_max, min_t_cnt=t_min, max_t_cnt=t_max, anchor=anchor, target=target)

print("Performing reconstruction...")
if args.core == 'CPU':
    reducer = umap.UMAP(n_components=2, 
                        n_neighbors = args.n_neighbors, 
                        min_dist = args.min_dist, 
                        n_epochs = args.n_epochs,
                        metric = args.metric,
                        
                        low_memory=False, 
                        verbose=True,
                        # local_connectivity = 30,
                        
                        random_state = 42)
    embedding = reducer.fit_transform(np.log1p(counts))

    # output reconstruction result
    a_recon = pd.DataFrame(embedding)
    a_recon.columns = ['xcoord','ycoord']
    a_recon.insert(loc=0, column=anchor, value=a_sel[anchor])
    a_recon.to_csv(os.path.join(out_dir, f'{anchor}_recon_loc.csv'), index=False)
        
    plot_shape(embedding, anchor, os.path.join(out_dir, f'UMAP_{anchor}.png'))
    plot_uniformness(embedding, counts, os.path.join(out_dir, f'UMAP_density_{anchor}.png'))
    plot_convex(embedding, anchor, os.path.join(out_dir, f'UMAP_convex_{anchor}.png'))

elif args.core == 'GPU':
    from cuml.manifold.umap import UMAP as cuUMAP
    min_dist_list = [0.3, 0.65, 0.99]
    n_epo_list = [10000, 50000, 100000]
    for md in min_dist_list:
        for nepo in n_epo_list:
            out_gpu = os.path.join(out_dir,'recon_GPU_{}_{}_{}'.format(m, md, nepo))
            if not os.path.exists(out_gpu):
                os.makedirs(out_gpu)
            
            reducer = cuUMAP(metric='cosine',
                                n_neighbors=25, 
                                min_dist=md,
                                n_components=2, 
                                verbose=True, 
                                n_epochs=nepo,
                                # local_connectivity = 30,
                                learning_rate = 1)
            embedding = reducer.fit_transform(np.log1p(counts))

            # output reconstruction result
            a_recon = pd.DataFrame(embedding)
            a_recon.columns = ['xcoord','ycoord']
            a_recon.insert(loc=0, column=anchor, value=a_sel[anchor])
            a_recon.to_csv(os.path.join(out_gpu, f'{anchor}_recon_loc.csv'), index=False)

            plot_shape(embedding, anchor, os.path.join(out_gpu, f'UMAP_{anchor}.png'))
            plot_uniformness(embedding, counts, os.path.join(out_gpu, f'UMAP_density_{anchor}.png'))
            plot_convex(embedding, anchor, os.path.join(out_gpu, f'UMAP_convex_{anchor}.png'))
else:
    assert False, f"ERROR: unknown core ({args.core})"

print("Done!")
