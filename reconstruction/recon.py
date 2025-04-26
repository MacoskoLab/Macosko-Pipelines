import os
import gc
import csv
import argparse
import numpy as np
import pandas as pd
import igraph as ig
import scipy.sparse as sp
from pypdf import PdfWriter
from helpers import *

def get_args():
    parser = argparse.ArgumentParser(description='2D embedding of reconstruction data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-c", "--cores", type=int, default=-1)
    parser.add_argument("-b", "--bead", type=int, default=2)
    parser.add_argument("-K", "--knn_filter", help="KNN graph filter", type=bool, default=False)
    
    parser.add_argument("-nn", "--n_neighbors", type=int, default=30)
    parser.add_argument("-md", "--min_dist", type=float, default=0.1)
    parser.add_argument("-s", "--spread", type=float, default=1.0)
    parser.add_argument("-lc", "--local_connectivity", type=int, default=1)
    parser.add_argument("-rs", "--repulsion_strength", type=float, default=1.0)
    parser.add_argument("-nsr", "--negative_sample_rate", type=int, default=10)
    parser.add_argument("-ne", "--n_epochs", type=int, default=500)
    
    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

args = get_args()

in_dir = args.in_dir             ; print(f"input directory: {in_dir}")
out_dir = args.out_dir           ; print(f"output directory: {out_dir}")
cores = args.cores               ; print(f"cores: {cores}")
bead = args.bead                 ; print(f"bead: {bead}")
knn_filter = args.knn_filter     ; print(f"KNN filter: {knn_filter}")

opts = dict()
opts["n_neighbors"] = args.n_neighbors
opts["min_dist"] = args.min_dist
opts["spread"] = args.spread
opts["local_connectivity"] = args.local_connectivity
opts["repulsion_strength"] = args.repulsion_strength
opts["negative_sample_rate"] = args.negative_sample_rate
opts["n_epochs"] = args.n_epochs
print("UMAP options:\n" + '\n'.join(f'    {k} = {v}' for k, v in opts.items()))
# learning_rate, set_op_mix_ratio, disconnection_distance

# Create programmatic output folder name
name = f"UMAP{bead}_"
name += "_".join([f"{''.join(s[0] for s in k.split('_'))}{v}" for k,v in opts.items()])
name += ("_K" if knn_filter else "")
print(f"name: {name}")

# Check in/out directories
assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['sb1.txt.gz', 'sb2.txt.gz', f'knn{bead}.npz'])
out_dir = os.path.join(out_dir, name)
if os.path.exists(out_dir):
    print(f"WARNING: output {out_dir} exists, overwriting...")
os.makedirs(out_dir, exist_ok=True)
assert os.path.exists(out_dir)

# Get number of allowed cores
avail_cores = len(os.sched_getaffinity(0))
if cores < 1:
    cores = avail_cores
if cores > avail_cores:
    print(f"WARNING: {cores} cores queried but only {avail_cores} available")
    cores = avail_cores
assert 0 < cores <= avail_cores
del avail_cores
print(f"Using {cores} cores")

### Pre-process KNN ############################################################

print("\nLoading KNN...")

# Read KNN matrix
knn_matrix = sp.load_npz(os.path.join(in_dir, f'knn{bead}.npz')).tocsr()
knn_matrix.setdiag(0) ; knn_matrix.eliminate_zeros()
assert np.diff(knn_matrix.indptr).max() >= opts["n_neighbors"] - 1

knn_mask = KNNMask(knn_matrix)

# Prune non-reciprocated edges
# m = np.abs(knn_matrix - knn_matrix.T) > np.min(knn_matrix.data)/2
# knn_matrix[m] = 0
# knn_matrix.eliminate_zeros()
# del m ; gc.collect()

print("\nFiltering beads with no neighbors...")
m = np.diff(knn_matrix.indptr) != 0
print(f"{np.sum(~m)} beads removed")
knn_matrix = knn_mask.apply_mask(knn_matrix, m, "No neighbors")
del m ; gc.collect()

print(f"\nReducing KNN to {opts['n_neighbors']} neighbors...")
knn_matrix = csr_k_nearest(knn_matrix, opts['n_neighbors'] - 1)
gc.collect()

print("\nRemoving disconnected components...")
n_components, labels = sp.csgraph.connected_components(knn_matrix, directed=False, return_labels=True)
m = labels == np.bincount(labels).argmax()
knn_matrix = knn_mask.apply_mask(knn_matrix, m, "Disconnected")
print(f"{np.sum(~m)} beads removed")
del m, n_components, labels ; gc.collect()

print("\nCreating KNN graph...")
G = ig.Graph(n = knn_matrix.shape[0],
             edges = zip(*knn_matrix.nonzero()),
             #edge_attrs = {'dist': knn_matrix.data, 'sim': 1-knn_matrix.data},
             directed = False)
G.simplify(multiple=True, loops=True, combine_edges="max")

print("\nPlotting KNN metrics...")
knn_indices, knn_dists = knn_matrix2indist(knn_matrix, opts["n_neighbors"])
nnd = knn_dists[:,1] ; nnd = nnd[np.isfinite(nnd)]
fnd = knn_dists[:,-1] ; fnd = fnd[np.isfinite(fnd)]
ies = np.array(G.degree(mode="in"))
tlu = np.array(G.transitivity_local_undirected(mode="zero", weights=None))
# con = np.array(G.constraint(weights="weight"))
# nss = np.array(G.neighborhood_size(order=2))
# sts = np.array(G.strength(weights="weight"))
fig = knn_plot(nnd, fnd, ies, tlu, opts['n_neighbors'])
fig.savefig(os.path.join(out_dir, "knn.pdf"), dpi=200)
del fig, knn_indices, knn_dists

if knn_filter:
    print("Running KNN filter...")

    # Filter low clustering coefficient (<0.2)
    m = tlu>=0.2
    knn_matrix = knn_mask.apply_mask(knn_matrix, m, "KNN filter")
    G = G.subgraph(np.where(m)[0])
    print(f"{np.sum(~m)} beads removed")
    del m ; gc.collect()
    
    # Filter beads with far nearest neighbors?
    # Filter too-high or too-low in-edges?
    # New collision detector?

# At this point we have G, knn_matrix, and knn_mask

### Initialization #############################################################

print("\nGenerating Leiden initialization...")

def leiden_init(G):
    import leidenalg as la

    print("    Leiden cluster...")
    # G.community_leiden(objective_function='modularity', weights=None, resolution=160, beta=0.01, n_iterations=2)
    # community_fastgreedy
    # community_infomap
    # community_label_propagation
    # community_leading_eigenvector
    # community_leading_eigenvector_naive
    # community_leiden
    # community_multilevel
    # community_optimal_modularity
    # community_spinglass
    # community_walktrap
    
    # https://leidenalg.readthedocs.io/en/stable/reference.html
    partition = la.find_partition(graph=G,
                                  partition_type=la.RBConfigurationVertexPartition,
                                  weights=None, # cosine?
                                  n_iterations=2,
                                  max_comm_size=1000,
                                  resolution_parameter=160)
    membership = np.array(partition.membership)
    print(f'    Number of clusters: {len(np.unique(membership))}')
    print(f'    Modularity: {partition.modularity}')

    # Contract KNN graph -> Leiden graph
    print("    Leiden graph...")
    G.contract_vertices(membership, combine_attrs=None)
    G.simplify(multiple=False, loops=True, combine_edges=None)
    G.es["weight"] = [1 - w for w in G.es["weight"]] # cosine distance -> similarity
    G.simplify(multiple=True, loops=True, combine_edges=sum)
    # Convert Leiden graph to CSR
    i, j = zip(*G.get_edgelist())
    ic_edges = sp.csr_matrix((G.es["weight"], (i,j)), shape=(G.vcount(), G.vcount()))
    ic_edges = ic_edges + ic_edges.T # ic_edges = ic_edges.maximum(ic_edges.T)
    ic_edges.setdiag(0)
    ic_edges.eliminate_zeros()

    # Create the Leiden KNN (for UMAP)
    print("    Leiden KNN...")
    def my_knn(mat):
        from sklearn.preprocessing import normalize
        from sparse_dot_topn import sp_matmul_topn
        mat_norm = normalize(mat, norm='l2', axis=1, copy=True)
        C = sp_matmul_topn(mat_norm, mat_norm.T, top_n=15, n_threads=1)
        C.data = 1 - C.data # cosine similarity -> distance
        C.data = np.maximum(C.data, 0)
        C.data = np.minimum(C.data, 1)
        C.setdiag(0)
        C.eliminate_zeros()
        # assert isinstance(C, sp.csr_matrix) and np.all(np.diff(C.indptr) == 14)
        return(C)
        
    knn_indices, knn_dists = knn_matrix2indist(my_knn(ic_edges))
    assert ic_edges.shape[0] == ic_edges.shape[1] == knn_indices.shape[0] == knn_dists.shape[0]

    print("    Leiden UMAP...")
    from umap import UMAP        
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   random_state = None,
                   verbose = False,
                   n_neighbors=15,
                   #min_dist=0.01,
                   n_epochs = 5000,
                   precomputed_knn = (knn_indices, knn_dists),
                   init = "spectral",
                   n_jobs = -1
                  )
    embedding = reducer.fit_transform(ic_edges)
    embedding -= np.mean(embedding, axis=0) # center
    return membership, ic_edges, embedding

if os.path.isfile(os.path.join(out_dir, 'init.npz')):
    # Loading previous
    print("TODO")
    membership, ic_edges, embedding = leiden_init(G)
else:
    membership, ic_edges, embedding = leiden_init(G)

print("    Leiden plots...")
# Metrics plot
# umis, connections, etc - leiden filter
# Initialization plot
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(embedding[:, 0], embedding[:, 1], s=1)
ax.set_title('Leiden Initialization')
ax.axis('equal')
fig.savefig(os.path.join(out_dir, "init.pdf"), dpi=200)
del fig, ax

print("    Leiden initialization...")
# Generate bead initialization - stacked
'''
init = embedding[membership]
'''
# Generate bead initialization - random jitter (default)
'''
init = embedding[membership]
from scipy.spatial import KDTree
sd = np.mean(KDTree(embedding).query(embedding, k=2)[0][:,1]) / 1000
init += np.random.normal(loc=0, scale=sd, size=(init.shape[0],2))
'''
# Generate bead initialization - (weighted?) average
knn_indices, knn_dists = knn_matrix2indist(knn_matrix, opts["n_neighbors"])
init = embedding[membership]
for _ in range(3):
    init = np.stack((np.mean(init[knn_indices, 0], axis=1),
                     np.mean(init[knn_indices, 1], axis=1)), axis=1)
del knn_indices, knn_dists

print("Saving Leiden initialization...")
np.savez_compressed(os.path.join(out_dir, 'init.npz'), embedding=embedding, membership=membership, init=init)
del G, knn_matrix ; gc.collect()

# at this point, we have:
#  init
#  final knn_mask (for sb)

### Load + subset data #########################################################

print(f"\nFiltering beads with <{opts['n_neighbors']} neighbors...")
while any(np.diff(knn_matrix.indptr) < opts['n_neighbors'] - 1):
    m = np.diff(knn_matrix.indptr) >= opts['n_neighbors'] - 1
    print(f"{np.sum(~m)} beads removed")
    knn_matrix = knn_mask.apply_mask(knn_matrix, m)
    del m ; gc.collect()
assert all(np.diff(knn_matrix.indptr) >= opts['n_neighbors'] - 1)

# Load barcodes
print('Loading the barcodes...')
sb = pd.read_csv(os.path.join(in_dir, f"sb{bead}.txt.gz"), header=None, names=['sb'])

# Subset barcodes
assert len(knn_mask.mask) == sb.shape[0]
assert type(sb) == pd.core.frame.DataFrame
sb = sb[knn_mask.mask]
gc.collect()

#assert mat.shape[0] == knn_indices.shape[0] == knn_dists.shape[0] == sb.shape[0]
assert knn_matrix.shape[0] == knn_matrix.shape[1] == sb.shape[0]

### UMAP TIME ##################################################################

print("\nRunning UMAP...")
embeddings = []
if opts["n_epochs"] <= 1000:
    embeddings.append(my_umap(knn_matrix, init, opts, cores))
else:
    n_epochs = opts["n_epochs"]
    opts["n_epochs"] = 1000
    embeddings.append(my_umap(knn_matrix, init, opts, cores))
    for i in range(int(np.ceil(n_epochs/1000))-1):
        embeddings.append(my_umap(knn_matrix, embeddings[-1], opts, cores))
    opts["n_epochs"] = n_epochs

### WRITE RESULTS ##############################################################

print("\nWriting results...")
embedding = embeddings[-1]

# Save the embeddings
np.savez_compressed(os.path.join(out_dir, "embeddings.npz"), *embeddings)

# Create the Puck file
assert embedding.shape[0] == len(sb)
pd.concat([sb, pd.DataFrame(embedding)], axis=1).to_csv(os.path.join(out_dir, 'Puck.csv'), index=False, header=False)

# Plot the umap
title = f"UMAP hexbin ({embedding.shape[0]:} beads) [{opts['n_epochs']} epochs]"
fig, ax = hexmap(embedding, title=title)
fig.savefig(os.path.join(out_dir, "umap.pdf"), dpi=200)

if len(embeddings) > 1:
    # Plot the intermediate embeddings
    fig, axes = hexmaps(embeddings, titles=[(i+1)*1000 for i in range(len(embeddings))])
    fig.savefig(os.path.join(out_dir, "umaps.pdf"), dpi=200)

    # Plot the convergence
    fig, axes = convergence_plot(embeddings)
    fig.savefig(os.path.join(out_dir, "convergence.pdf"), dpi=200) ; del fig

# Plot the weighted embeddings
title = f"UMAP hexbin (logUMI-weighted)"
weights = np.log1p(mat).sum(axis=1).A1
fig, ax = hexmap(embedding, weights, title=title)
fig.savefig(os.path.join(out_dir, "umap_umi.pdf"), dpi=200)

# Make summary pdf
# creation order: knn_filter, knn_graph, knn, leiden, umap, umaps, convergence, umap_umi
names = ["umap", "convergence", "umaps", "umap_umi", "init", "knn"]
paths = [os.path.join(out_dir, n+".pdf") for n in names]
files = [p for p in paths if os.path.isfile(p)]
if len(files) > 0:
    merger = PdfWriter()

    for file_name in files:
        merger.append(file_name)

    merger.write(os.path.join(out_dir, "summary.pdf"))
    merger.close()
    [os.remove(file) for file in files]

print("Done!")
