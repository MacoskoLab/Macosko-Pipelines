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
# os.chdir("/home/nsachdev/reconstruction/1.2")
assert all(os.path.isfile(file) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz'])

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-n", "--n_neighbors", help="the number of neighboring sample points used for manifold approximation", type=int, default=25)
    parser.add_argument("-m", "--metric", help="the metric to use to compute distances in high dimensional space", type=str, default="cosine")
    parser.add_argument("-d", "--min_dist", help="the effective minimum distance between embedded points", type=float, default=0.99)
    parser.add_argument("-I", "--init", help="how to initialize the low dimensional embedding", type=str, default="spectral")
    parser.add_argument("-N", "--n_epochs", help="the number of training epochs to be used in optimizing the low dimensional embedding", type=int, default=20000)
    # parser.add_argument("-c", "--core", help="define core type to use (CPU or GPU)", type=str, default="CPU")
    args = parser.parse_args()
    return args

args = get_args()
n_neighbors = args.n_neighbors ; print(f"n_neighbors = {n_neighbors}")
metric = args.metric           ; print(f"metric = {metric}")
min_dist = args.min_dist       ; print(f"min_dist = {min_dist}")
init = args.init               ; print(f"init = {init}")
n_epochs = args.n_epochs       ; print(f"n_epochs = {n_epochs}")
out_dir = f"ANCHOR_{n_neighbors}_{metric}_{min_dist}_{init}_{n_epochs}"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("\nReading the matrix...")
df = pd.read_csv('matrix.csv.gz', compression='gzip', header=None, names=['sb1', 'sb2', 'umi'])
df.sb1 -= 1 # convert from 1- to 0-indexed
df.sb2 -= 1 # convert from 1- to 0-indexed
with gzip.open('sb1.txt.gz', 'rt') as f:
    sb1 = [line.strip() for line in f.readlines()]
with gzip.open('sb2.txt.gz', 'rt') as f:
    sb2 = [line.strip() for line in f.readlines()]
assert sorted(list(set(df.sb1))) == list(range(len(sb1)))
assert sorted(list(set(df.sb2))) == list(range(len(sb2)))
print(f"{len(sb1)} R1 barcodes")
print(f"{len(sb2)} R2 barcodes")

# ???
grouped_sum = df.groupby('sb1')['umi'].sum()
sb1_to_keep = grouped_sum[grouped_sum <= 7000].index
df = df[df['sb1'].isin(sb1_to_keep)]

# Rows are the anchor beads I wish to recon
# Columns are the features used for judging similarity
mat = coo_matrix((df['umi'], (df['sb2'], df['sb1'])))

print("\nComputing the KNN...")
knn = nearest_neighbors(mat,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        metric_kwds=None,
                        angular=False,
                        random_state=None,
                        low_memory=False,
                        verbose=True
                       )

def my_umap(mat, n_epochs, init=init):
    reducer = UMAP(n_components = 2,
               random_state=None,
               low_memory = False,
               verbose = True,
               precomputed_knn = knn,
               metric = metric,
               n_neighbors = n_neighbors,
               init = init,
               min_dist = min_dist,
               n_epochs = n_epochs,
               )
    embedding = reducer.fit_transform(np.log1p(mat))
    return(embedding)

print("\nRunning UMAP...")
embeddings = []
embeddings.append(my_umap(mat, n_epochs=10))
embeddings.append(my_umap(mat, n_epochs=100, init=embeddings[-1]))
for i in range(round(n_epochs/1000)):
    print(i)
    embeddings.append(my_umap(mat, init=embeddings[-1], n_epochs=1000))

print("\ndone")
np.savez(os.path.join(out_dir,"embeddings.npz"), *embeddings)
