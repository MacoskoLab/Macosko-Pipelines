import os
import gc
import sys
import gzip
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from umap import UMAP
import umap.plot
from umap.umap_ import nearest_neighbors

# os.chdir("/home/nsachdev/recon/data/615-2cm")
# os.chdir("/home/nsachdev/recon/data/609-6mm")

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    # parser.add_argument("-c", "--core", help="define core type to use (CPU or GPU)", type=str, default="CPU")
    # bead type? tags or seq? "-e", "--exptype", help="define experiment type (seq or tags)", type=str, required=True,

    parser.add_argument("-a", "--algorithm", help="dimensionality reduction algo", type=str, default="umap")
    
    parser.add_argument("-n", "--n_neighbors", help="the number of neighboring sample points used for manifold approximation", type=int, default=25)
    parser.add_argument("-d", "--min_dist", help="the effective minimum distance between embedded points", type=float, default=0.99)
    parser.add_argument("-s", "--spread", help="the effective scale of embedded points", type=float, default=1.0)
    parser.add_argument("-N", "--n_epochs", help="the number of training epochs to be used in optimizing the low dimensional embedding", type=int, default=10000)
    parser.add_argument("-I", "--init", help="how to initialize the low dimensional embedding", type=str, default="spectral")
    parser.add_argument("-m", "--metric", help="the metric to use to compute distances in high dimensional space", type=str, default="cosine")
    
    args, unknown = parser.parse_known_args()
    return args

args = get_args()
algo = args.algorithm ; print(f"algorithm = {algo}")
if algo == "umap":
    n_neighbors = args.n_neighbors ; print(f"n_neighbors = {n_neighbors}")
    min_dist = args.min_dist       ; print(f"min_dist = {min_dist}")
    spread = args.spread           ; print(f"spread = {spread}")
    n_epochs = args.n_epochs       ; print(f"n_epochs = {n_epochs}")
    init = args.init               ; print(f"init = {init}")
    metric = args.metric           ; print(f"metric = {metric}")
    out_dir = os.path.join(args.out_dir, f"ANCHOR_n={n_neighbors}_d={min_dist}_s={spread}_N={n_epochs}_I={init}_m={metric}")
else:
    out_dir = os.path.join(args.out_dir, f"ANCHOR")

in_dir = args.in_dir
assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz'])

if not os.path.exists(out_dir):
        os.makedirs(out_dir)

print("\nReading the matrix...")
df = pd.read_csv(os.path.join(in_dir, 'matrix.csv.gz'), compression='gzip', header=None, names=['sb1', 'sb2', 'umi'])
df.sb1 -= 1 # convert from 1- to 0-indexed
df.sb2 -= 1 # convert from 1- to 0-indexed
with gzip.open(os.path.join(in_dir, 'sb1.txt.gz'), 'rt') as f:
    sb1 = [line.strip() for line in f.readlines()]
with gzip.open(os.path.join(in_dir, 'sb2.txt.gz'), 'rt') as f:
    sb2 = [line.strip() for line in f.readlines()]
assert sorted(list(set(df.sb1))) == list(range(len(sb1)))
assert sorted(list(set(df.sb2))) == list(range(len(sb2)))
print(f"{len(sb1)} R1 barcodes")
print(f"{len(sb2)} R2 barcodes")

# Rows are the anchor beads I wish to recon
# Columns are the features used for judging similarity
mat = coo_matrix((df['umi'], (df['sb2'], df['sb1'])))

### PLOTS ######################################################################

# fig, axs = plt.subplots(3, 2, figsize=(7, 8))

# def plot_histogram(ax, df, col, operation):
#     if operation == 'umi':
#         gdf = df.groupby(col)['umi'].sum().reset_index()
#         gdf["umi"] = np.log10(gdf["umi"])
#         xlab = "log10 umi"
#     elif operation == 'connections':
#         gdf = df.groupby(col)['umi'].size().reset_index()
#         gdf.columns = [col, 'umi']
#         gdf["umi"] = np.log10(gdf["umi"])
#         xlab = "log10 connections"
#     elif operation == 'entropy':
#         gdf = df.groupby(col)['umi'].apply(shannon_entropy).reset_index()
#         gdf.columns = [col, 'umi']
#         xlab = "entropy"

#     #meanumi = np.log10(gdf['umi'].mean())
#     #meanumi = gdf['umi'].mean()
#     ax.hist(gdf['umi'], bins=100)
#     #ax.axvline(meanumi, color='red', linestyle='dashed', linewidth=2)
#     #ax.text(meanumi, ax.get_ylim()[1]*0.9, f'Mean: {10**meanumi:.1f}', color='red', fontsize=12, ha='center')
#     ax.set_yscale('log')
#     ax.set_title(f"Histogram of total {operation} per {col}")
#     ax.set_xlabel(xlab)
#     ax.set_ylabel('log10 count')

# plot_histogram(axs[0,0], df, 'sb1', 'umi')
# plot_histogram(axs[1,0], df, 'sb1', 'connections')
# #plot_histogram(axs[2,0], df, 'sb1', 'entropy')
# plot_histogram(axs[0,1], df, 'sb2', 'umi')
# plot_histogram(axs[1,1], df, 'sb2', 'connections')
# #plot_histogram(axs[2,1], df, 'sb2', 'entropy')

# # Adjust layout to prevent overlap
# plt.tight_layout()

# # Save the plot as a PDF
# plt.savefig('umi_histograms.pdf')

### UMAP TIME ##################################################################

def my_umap(mat, n_epochs, init=init):
    reducer = UMAP(n_components = 2,
                   random_state = None,
                   low_memory = True,
                   verbose = True,
                   precomputed_knn = knn,
                   
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   spread = spread,
                   n_epochs = n_epochs,
                   init = init,
                   metric = metric
                  )
    embedding = reducer.fit_transform(np.log1p(mat))
    return(embedding)

if algo == "umap":
    
    print("\nComputing the KNN...")
    knn = nearest_neighbors(mat,
                            n_neighbors=n_neighbors,
                            metric=metric,
                            metric_kwds=None,
                            angular=False,
                            random_state=None,
                            low_memory=True,
                            use_pynndescent=True,
                            n_jobs=-1,
                            verbose=True
                           )

    print("\nRunning UMAP...")
    embeddings = []
    embeddings.append(my_umap(mat, n_epochs=10))
    embeddings.append(my_umap(mat, n_epochs=100, init=embeddings[-1]))
    embeddings.append(my_umap(mat, n_epochs=890, init=embeddings[-1]))
    for i in range(round(n_epochs/1000)):
        print(i)
        embeddings.append(my_umap(mat, init=embeddings[-1], n_epochs=1000))
    
    print("\ndone")
    np.savez(os.path.join(out_dir, "embeddings.npz"), *embeddings)
