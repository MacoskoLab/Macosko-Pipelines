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

# Load data ####################################################################

print("\nLoading data...")

# Load KNN matrix
knn_matrix = sp.load_npz(os.path.join(in_dir, f'knn{bead}.npz')).tocsr()
knn_matrix.setdiag(0) ; knn_matrix.eliminate_zeros()
assert np.diff(knn_matrix.indptr).max() >= opts["n_neighbors"] - 1

# Load sb
sb = pd.read_csv(os.path.join(in_dir, f"sb{bead}.txt.gz"), header=None, names=['sb'])['sb'].to_numpy()
assert len(sb) == knn_matrix.shape[0] == knn_matrix.shape[1]

# Load the number of UMIs per bead
if os.path.isfile(os.path.join(in_dir, f"readumi_per_sb{bead}.csv.gz")):
    df = pd.read_csv(os.path.join(in_dir, f"readumi_per_sb{bead}.csv.gz"),
                     sep=',', compression='gzip', usecols=[f'sb{bead}', 'umis'])
    df = df.rename(columns={f'sb{bead}': 'sb'})
    umi = df.set_index('sb').loc[sb, "umis"].to_numpy()
    del df
else:
    df = pd.read_csv(os.path.join(in_dir, f"matrix.csv.gz"),
                     sep=',', compression='gzip', usecols=[f'sb{bead}_index', 'umi'])
    df = df.rename(columns={f'sb{bead}_index': 'sb_index'})
    df.sb_index -= 1
    df = df.groupby('sb_index', as_index=False)['umi'].sum()
    df = df.sort_values('sb_index')
    umi = df['umi'].to_numpy()
    del df

knn = KNN(knn_matrix, sb, umi)
del knn_matrix, sb, umi ; gc.collect()

print("\nFiltering beads with no neighbors...")
knn.apply_mask(np.diff(knn.matrix.indptr) != 0, "No neighbors")

### Load KNN graph ############################################################

print("\nCreating KNN graph...")
knn_matrix_k = csr_k_nearest(knn.matrix, k=opts['n_neighbors']-1)
G = ig.Graph(n = knn_matrix_k.shape[0],
             edges = zip(*knn_matrix_k.nonzero()),
             edge_attrs = {'weight': 1-knn_matrix_k.data}, # cosine similarity weights
             directed = False)
G.simplify(multiple=True, loops=True, combine_edges="max")
del knn_matrix_k ; gc.collect()

print("\nPlotting KNN metrics...")
knn_indices, knn_dists = knn_matrix2indist(knn.matrix, k=opts["n_neighbors"])
nnd = knn_dists[:,1] ; nnd = nnd[np.isfinite(nnd)]
fnd = knn_dists[:,-1] ; fnd = fnd[np.isfinite(fnd)]
ies = np.array(G.degree(mode="in"))
tlu = np.array(G.transitivity_local_undirected(mode="zero", weights='weight'))
sts = np.array(G.strength(weights="weight"))
con = np.array(G.constraint(weights="weight"))

def my_histogram(ax, vec, title, xlab):
    ax.hist(vec, bins=100)
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel("Beads")
    
    xmean = np.mean(vec)
    ax.axvline(xmean, color='red', linestyle='dotted', linewidth=1.5)
    ax.text(xmean + ax.get_xlim()[1] * 0.02, ax.get_ylim()[1] * 0.925,
            f'Mean: {xmean:.2f}', color='red', ha='left')

fig, axs = plt.subplots(3, 2, figsize=(8, 8))
my_histogram(axs[0,0], nnd, "Nearest neighbor distance", "Cosine distance")
my_histogram(axs[0,1], fnd, f"Furthest neighbor ({opts['n_neighbors']}) distance", "Cosine distance")
my_histogram(axs[1,0], ies, f"Number of in-edges", "In-edges")
my_histogram(axs[1,1], tlu, f"Clustering coefficient", "transitivity_local_undirected")
my_histogram(axs[2,0], sts, f"Strength", "strength")
my_histogram(axs[2,1], con, f"Burt's constraint", "constraint")
plt.tight_layout()
fig.savefig(os.path.join(out_dir, "knn.pdf"), dpi=200)
del fig, axs, knn_indices, knn_dists, my_histogram

if knn_filter:
    print("Running KNN filter...")

    # Filter low clustering coefficient (<0.2)
    m = tlu>=0.2
    knn.apply_mask(m, "KNN filter")
    G = G.subgraph(np.where(m)[0])
    del m ; gc.collect()
    
    # Filter beads with far nearest neighbors?
    # Filter too-high or too-low in-edges?
    # New collision detector?

### Compute Initialization #####################################################

def leiden_init(G):
    import leidenalg as la

    print("    Leiden cluster...")
    # https://igraph.org/python/api/0.9.7/igraph._igraph.GraphBase.html#community_edge_betweenness
    # G.community_leiden(objective_function='modularity', weights=None, resolution=160, beta=0.01, n_iterations=2)
    # G.community_edge_betweenness()
    # G.community_fastgreedy
    # G.community_infomap
    # G.community_label_propagation
    # G.community_leading_eigenvector
    # G.community_leading_eigenvector_naive
    # G.community_leiden
    # G.community_multilevel
    # G.community_optimal_modularity
    # G.community_spinglass
    # G.community_walktrap
    
    # https://leidenalg.readthedocs.io/en/stable/reference.html
    partition = la.find_partition(graph=G,
                                  partition_type=la.RBConfigurationVertexPartition,
                                  weights='weight',
                                  n_iterations=2,
                                  #max_comm_size=1000,
                                  resolution_parameter=160)
    membership = np.array(partition.membership)
    print(f'    Number of clusters: {len(np.unique(membership))}')
    print(f'    Modularity: {partition.modularity}')

    # Contract KNN graph -> Leiden graph
    print("    Leiden graph...")
    G.contract_vertices(membership, combine_attrs=None)
    G.simplify(multiple=True, loops=True, combine_edges=sum)
    
    # Convert Leiden graph to symmetric CSR
    i, j = zip(*G.get_edgelist())
    ic_edges = sp.csr_matrix((G.es["weight"], (i,j)), shape=(G.vcount(), G.vcount()))
    ic_edges = ic_edges.maximum(ic_edges.T)
    assert ic_edges.has_canonical_format
    assert np.all(ic_edges.diagonal() == 0)
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

print("\nGenerating Leiden embedding...")
assert G.vcount() == knn_matrix.shape[0] == knn_matrix.shape[1]
membership, ic_edges, embedding = leiden_init(G)
assert membership.shape[0] == knn_matrix.shape[0] == knn_matrix.shape[1]
assert G.vcount() == embedding.shape[0] == np.unique(membership).shape[0]
del G ; gc.collect()

### Process Initialization #####################################################

print("\nLeiden plots...")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(x=embedding[:, 0], y=embedding[:, 1], s=1)
ax.set_title('Leiden Initialization')
ax.axis('equal')
fig.savefig(os.path.join(out_dir, "init.pdf"), dpi=200)
del fig, ax

fig, axs = plt.subplots(2, 2, figsize=(8, 8))
def my_plot(embedding, vec, lab, axs):
    axs[0].hist(vec, bins=30, color='skyblue', edgecolor='black')
    axs[0].set_title(lab)

    axs[1].scatter(x=embedding[:,0], y=embedding[:,1], c=vec, cmap='viridis', s=2)
    axs[1].set(xticks=[], yticks=[], xticklabels=[], yticklabels=[])
    [spine.set_visible(False) for spine in axs[1].spines.values()]
    axs[1].set_aspect('equal')
    axs[1].set_title(lab)

my_plot(embedding, np.bincount(membership), '#Beads', axs[0])
my_plot(embedding, np.sum(ic_edges, axis=0).A1, '#Edges', axs[1])
fig.savefig(os.path.join(out_dir, "init2.pdf"), dpi=200)
del fig, axs, my_plot, ic_edges

print("\nLeiden filter...")
# TODO
# (make sure to update knn_mask/knn_matrix and membership/ic_edges/embedding)

print(f"\nFiltering beads with <{opts['n_neighbors']} neighbors...")
while any(np.diff(knn_matrix.indptr) < opts['n_neighbors'] - 1):
    m = np.diff(knn_matrix.indptr) >= opts['n_neighbors'] - 1
    knn_matrix = knn_mask.apply_mask(knn_matrix, m, "Too few neighbors")
    membership = membership[m]
    print(f"{np.sum(~m)} beads removed")
    del m ; gc.collect()
assert all(np.diff(knn_matrix.indptr) >= opts['n_neighbors'] - 1)
assert knn_matrix.shape[0] == knn_matrix.shape[1] == membership.shape[0]

print("\nRemoving disconnected components...")
knn_matrix_k = csr_k_nearest(knn_matrix, k=opts['n_neighbors']-1)
n_components, labels = sp.csgraph.connected_components(knn_matrix_k, directed=False, return_labels=True)
m = labels == np.bincount(labels).argmax()
knn_matrix = knn_mask.apply_mask(knn_matrix, m, "Disconnected")
membership = membership[m]
print(f"{np.sum(~m)} beads removed")
del m, n_components, labels, knn_matrix_k ; gc.collect()
assert knn_matrix.shape[0] == knn_matrix.shape[1] == membership.shape[0]

### NO MORE FILTERING PAST THIS POINT ###

print("\nLeiden initialization...")
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
# Generate bead initialization - weighted average
knn_indices, knn_dists = knn_matrix2indist(knn_matrix, opts["n_neighbors"])
weights = 1-knn_dists ; weights[knn_indices < 0] = 0
init = embedding[membership]
for _ in range(3):
    init = np.stack((np.average(init[knn_indices, 0], axis=1, weights=weights),
                     np.average(init[knn_indices, 1], axis=1, weights=weights)), axis=1)
del knn_indices, knn_dists, weights

print("Saving Leiden initialization...")
np.savez_compressed(os.path.join(out_dir, 'init.npz'), embedding=embedding, membership=membership, init=init)
assert np.unique(membership).shape[0] <= embedding.shape[0] and np.max(membership) < embedding.shape[0]
assert init.shape[0] == knn_matrix.shape[0] == knn_matrix.shape[1]
del membership, embedding ; gc.collect()

### UMAP TIME ##################################################################

print("\nRunning UMAP...")
from umap import UMAP
reducer = UMAP(n_components = 2,
               metric = "precomputed",
               random_state = None,
               verbose = True,
               low_memory=True,
               init = init,
               n_jobs = cores,
               **opts
              )
embedding = reducer.fit_transform(knn_matrix.maximum(knn_matrix.T))
assert np.all(np.isfinite(embedding))

### WRITE RESULTS ##############################################################

print("\nWriting results...")

# Create the Puck file
sb = sb[knn_mask.mask].reset_index(drop=True)
assert embedding.shape[0] == len(sb)
pd.concat([sb, pd.DataFrame(embedding)], axis=1).to_csv(os.path.join(out_dir, 'Puck.csv'), index=False, header=False)

# Plot the umap
title = f"UMAP hexbin ({embedding.shape[0]:} beads) [{opts['n_epochs']} epochs]"
fig, ax = hexmap(embedding, title=title)
fig.savefig(os.path.join(out_dir, "umap.pdf"), dpi=200)

# Plot the weighted embeddings
# title = f"UMAP hexbin (logUMI-weighted)"
# weights = np.log1p(mat).sum(axis=1).A1
# fig, ax = hexmap(embedding, weights, title=title)
# fig.savefig(os.path.join(out_dir, "umap_umi.pdf"), dpi=200)

# Make summary pdf
# creation order: knn,  init, init2, umap, umap_umi
names = ["umap", "umap_umi", "knn", "init", "init2"]
paths = [os.path.join(out_dir, n+".pdf") for n in names]
files = [p for p in paths if os.path.isfile(p)]
if len(files) > 0:
    merger = PdfWriter()

    for file_name in files:
        merger.append(file_name)

    merger.write(os.path.join(out_dir, "summary.pdf"))
    merger.close()
    #[os.remove(file) for file in files]

print("Done!")
