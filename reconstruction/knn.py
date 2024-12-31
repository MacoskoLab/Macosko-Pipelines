import os
import gc
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.preprocessing import normalize
from sparse_dot_topn import sp_matmul_topn, zip_sp_matmul_topn

def get_args():
    parser = argparse.ArgumentParser(description='Run KNN on a diffusion matrix')
    parser.add_argument("-i", "--in_dir", type=str, default=".")
    parser.add_argument("-o", "--out_dir", type=str, default=".")
    parser.add_argument("-n", "--n_neighbors", type=int, default=150)
    parser.add_argument("-b", "--bead", type=int, default=2)
    parser.add_argument("-c", "--cores", type=int, default=-1)
    parser.add_argument("-k", "--chunks", type=int, default=2)
    
    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

# Load arguments
args = get_args()
in_dir = args.in_dir           ; print(f"input directory = {in_dir}")
out_dir = args.out_dir         ; print(f"output directory = {out_dir}")
n_neighbors = args.n_neighbors ; print(f"n_neighbors = {n_neighbors}")
bead = args.bead               ; print(f"bead = {bead}")
cores = args.cores             ; print(f"cores = {cores}")
chunks = args.chunks           ; print(f"chunks = {chunks}")

# Check in/out directories
assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz'])
os.makedirs(out_dir, exist_ok=True)
assert os.path.exists(out_dir)

# Get number of allowed cores
avail_cores = len(os.sched_getaffinity(0))
if cores < 1:
    cores = avail_cores
if cores > avail_cores:
    print(f"WARNING: {cores} cores queried but only {avail_cores} available")
    cores = avail_cores

# Check remaining arguments
assert n_neighbors > 0
assert bead in [1, 2]
assert 0 < cores <= avail_cores
assert chunks > 0

# Read matrix file
print('Loading the matrix...')
df = pd.read_csv(os.path.join(in_dir, 'matrix.csv.gz'), compression='gzip')
df.sb1_index -= 1 # convert from 1- to 0-indexed
df.sb2_index -= 1 # convert from 1- to 0-indexed

if bead == 1:
    mat = sp.coo_matrix((df['umi'], (df['sb1_index'], df['sb2_index']))).tocsr()
    out_file = os.path.join(out_dir, f"knn1.npz")
if bead == 2:
    mat = sp.coo_matrix((df['umi'], (df['sb2_index'], df['sb1_index']))).tocsr()
    out_file = os.path.join(out_dir, f"knn2.npz")

if os.path.exists(out_file):
    print(f"WARNING: file {out_file} exists, overwriting...")

del df ; gc.collect()

# Pre-process the matrix
print('Pre-processing the matrix...')
def process(mat):
    mat_norm = np.log1p(mat.copy())
    mat_norm = normalize(mat_norm, norm='l2', axis=1, copy=True)
    mat_norm = mat_norm.astype(np.float32, copy=True)
    mat_norm.eliminate_zeros()
    return(mat_norm)

mat_norm = process(mat)
del mat ; gc.collect()

# Chunk the matrix
print('Chunking the matrix...')
ABs = [mat_norm[chunk] for chunk in np.array_split(range(mat_norm.shape[0]), chunks)]
del mat_norm ; gc.collect()
for AB in ABs:
    assert AB.getnnz() < 2**31-1, "Matrix is too large, increase 'chunks'"

# Calculate Top-N KNN
print('Multiplying sub-matrices...')
Cs = [[sp_matmul_topn(A, B.T,
                      top_n=n_neighbors,
                      n_threads=cores
                     ) for B in ABs] for A in ABs]
del ABs ; gc.collect()
# TODO: investigate 'sp_matmul_topn' parameters 'threshold' and 'density'

print('Combining sub-matrices...')
C = sp.vstack([zip_sp_matmul_topn(top_n=n_neighbors, C_mats=Cis) for Cis in Cs],
              format="csr", dtype=np.float32)
del Cs ; gc.collect()

print('Finalizing result...')
C.data = 1 - C.data
C.data = np.maximum(C.data, 0)
C.data = np.minimum(C.data, 1)
C.setdiag(0)
C.eliminate_zeros()

print('Saving output...')
sp.save_npz(out_file, C, compressed=True)

print("Done!")

### Archive: NNDescent #########################################################

# # Compute the KNN using NNDescent
# def knn_descent(mat, n_neighbors, metric="cosine", n_jobs=-1):    
#     from pynndescent import NNDescent
#     knn_search_index = NNDescent(
#         data=mat,
#         metric=metric,
#         metric_kwds={},
#         n_neighbors=n_neighbors,
#         n_trees=64, # originally None
#         # leaf_size=None,
#         pruning_degree_multiplier=3.0, # originally 1.5
#         diversify_prob=0.0, # originally 1.0
#         # tree_init=True,
#         # init_graph=None,
#         # init_dist=None,
#         random_state=None,
#         low_memory=True, # originally False
#         max_candidates=60, # originally None
#         max_rptree_depth=999999, # originally 200
#         n_iters=512, # originally None
#         delta=0.0001, # originally 0.001
#         n_jobs=n_jobs,
#         # compressed=False,
#         # parallel_batch_queries=False,
#         verbose=True # originally False
#     )
#     knn_indices, knn_dists = knn_search_index.neighbor_graph
#     return knn_indices, knn_dists

# # UMAP NNDescent wrapper
# def nearest_neighbors(mat, n_neighbors, metric="cosine", n_jobs=-1):
#     from umap.umap_ import nearest_neighbors
#     knn_indices, knn_dists, _ = nearest_neighbors(mat,
#                                     n_neighbors = n_neighbors,
#                                     metric = metric,
#                                     metric_kwds = {},
#                                     angular = False, # Does nothing?
#                                     random_state = None,
#                                     low_memory = True, # False?
#                                     use_pynndescent = True, # Does nothing?
#                                     n_jobs = n_jobs,
#                                     verbose = True
#                                 )
#     return knn_indices, knn_dists

# # Take the best neighbors between two NNDescent runs
# # knn_indices1 is row-sliced, knn_indices2 is original
# def knn_merge(knn_indices1, knn_dists1, knn_indices2, knn_dists2):
#     assert knn_indices1.shape == knn_dists1.shape
#     assert knn_indices2.shape == knn_dists2.shape
#     assert knn_indices1.shape[0] == knn_dists1.shape[0] == knn_indices2.shape[0] == knn_dists2.shape[0]
#
#     assert all(knn_indices2[:,0] == np.arange(knn_indices2.shape[0]))
#     assert np.all(knn_indices2 >= 0) and np.all(knn_dists2 >= 0)
#     index_map = dict(zip(knn_indices1[:,0], knn_indices2[:,0]))
#     knn_indices1 = np.vectorize(lambda i: index_map.get(i,-1))(knn_indices1)
#     assert all(knn_indices1[:,0] == knn_indices2[:,0])
#
#     k = max(knn_indices1.shape[1], knn_indices2.shape[1])
#     knn_indices = np.zeros((knn_indices1.shape[0], k), dtype=np.int32)
#     knn_dists = np.zeros((knn_indices1.shape[0], k), dtype=np.float64)
#
#     for i in range(knn_indices1.shape[0]):
#         d1 = {i:d for i,d in zip(knn_indices1[i], knn_dists1[i]) if i >= 0}
#         d2 = {i:d for i,d in zip(knn_indices2[i], knn_dists2[i]) if i >= 0}
#         d = d1 | d2
#
#         row = sorted(d.items(), key=lambda item: item[1])[:k]
#         inds, dists = zip(*row)
#
#         knn_indices[i,:] = inds
#         knn_dists[i,:] = dists
#
#     return knn_indices, knn_dists

# knn_indices, knn_dists = knn_descent(np.log1p(mat), n_neighbors, metric="cosine", n_jobs=cores)
