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

    parser.add_argument("-K", "--knn_filter", help="k-nearest neighbor filter", type=bool, default=False)
    parser.add_argument("-M", "--mnn", help="mutual nearest neighbors", type=bool, default=False)
    
    parser.add_argument("-nn", "--n_neighbors", type=int, default=45)
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
mnn = args.mnn                   ; print(f"MNN: {mnn}")

opts = dict()
opts["n_neighbors"] = args.n_neighbors
opts["min_dist"] = args.min_dist
opts["spread"] = args.spread
opts["local_connectivity"] = args.local_connectivity
opts["repulsion_strength"] = args.repulsion_strength
opts["negative_sample_rate"] = args.negative_sample_rate
opts["n_epochs"] = args.n_epochs
print(f"UMAP options:\n{'\n'.join(f'    {k} = {v}' for k, v in opts.items())}")
# learning_rate, set_op_mix_ratio, disconnection_distance

# Create programmatic output folder name
name = f"UMAP_sb{bead}_"
name += "_".join([f"{''.join(s[0] for s in k.split('_'))}{v}" for k,v in opts.items()])
name += "_" if mnn or knn_filter else ""
name += ("M" if mnn else "") + ("K" if knn_filter else "")
print(f"name: {name}")

# Check in/out directories
assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz', f'knn{bead}.npz'])
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
k = np.diff(knn_matrix.indptr).max()+1 # +1 because K includes self-node
assert opts["n_neighbors"] <= k

knn_mask = KNNMask(knn_matrix)

# Prune non-reciprocated edges
if mnn:
    print(f"Running MNN with {k} neighbors...")
    m = np.abs(knn_matrix - knn_matrix.T) > np.min(knn_matrix.data)/2
    knn_matrix[m] = 0
    knn_matrix.eliminate_zeros()
    del m ; gc.collect()

print("Filtering beads with no neighbors...")
m = np.diff(knn_matrix.indptr) != 0
print(f"{sum(~m)} beads removed")
knn_matrix = knn_mask.apply_mask(knn_matrix, m)
del m ; gc.collect()

print(f"Reducing KNN from {k} to {opts['n_neighbors']} neighbors...")
knn_matrix = csr_k_nearest(knn_matrix, opts['n_neighbors']-1)
gc.collect()

print("Creating KNN graph...")
G = ig.Graph(n = knn_matrix.shape[0],
             edges = zip(*knn_matrix.nonzero()),
             edge_attrs = {'weight': knn_matrix.data},
             directed = False)
G.simplify(combine_edges="max")

if knn_filter:
    print("Running KNN filter...")
    assert False, "TODO: incorporate new collision detector"
    #cc = G.connected_components()
    #tlu = G.transitivity_local_undirected(mode="zero", weights="weight")

    # # Remove weakly-connected components
    # print("Removing WCCs...")
    # n_components, labels = sp.csgraph.connected_components(knn_matrix, directed=True, connection='strong', return_labels=True)
    # m = labels == np.bincount(labels).argmax()
    # knn_matrix = knn_mask.apply_mask(knn_matrix, m)
    # print(f"{np.sum(~m)} WCC beads removed")
    # del m ; gc.collect()

    # # Filter low clustering coefficient
    # clustering = nx.clustering(G, nodes=None, weight=None)
    # data = [clustering[key] for key in sorted(clustering.keys())]
    # hist_z(axes[1,1], data, z_low=z_low)
    # axes[1,1].set_xlabel('Clustering coefficient')
    # axes[1,1].set_title('Local clustering coefficient')

    # Filter beads with far nearest neighbors (?)
    # Filter too-high or too-low in-edges (?)
    # this method needs to update knn_mask, knn_matrix, and G
    # makes sure the node order in G continues to match knn_matrix
    # print("X total beads removed")

knn_indices, knn_dists = knn_matrix2indist(knn_matrix)
del knn_matrix ; gc.collect()

# KNN graph plot (TODO)
tlu = G.transitivity_local_undirected(mode="zero", weights="weight")

# KNN indist plot
fig, ax = knn_plot(knn_indices, knn_dists, tlu)
fig.savefig(os.path.join(out_dir, "knn.pdf"), dpi=200)
del fig, ax

# Unused metrics
# G.strength(G.vs, weights="weight")
# G.reciprocity()
# G.personalized_pagerank()
# G.path_length_hist()
# G.neighborhood_size()
# G.neighborhood()
# G.knn()
# G.hub_score()
# G.diversity()
# G.degree()
# G.coreness()
# G.constraint()
# G.cocitation()
# G.closeness()
# G.cliques()
# G.bibcoupling()
# G.betweenness()
# G.authority_score()
# G.assortativity()

# at this point, we have:
#   final G (for leiden)
#   final knn_indices/knn_dists (for umap)
#   final knn_mask (for matrix, sb2)

### Initialization #############################################################

print("\nGenerating Leiden initialization...")

def leiden_init(G):
    import leidenalg as la

    print("    Leiden cluster...")
    # https://igraph.org/python/doc/api/igraph.Graph.html#community_leiden
    #a1 = G.community_leiden(objective_function='modularity', weights=None, resolution=160, beta=0.01, n_iterations=2, node_weights=None)
    # CPM, res=1 - 11.39
    # mod, res=1 - 33.968
    # CPM, res=160 - 11.588
    # mod, res=160 - 32.728
    # start_time = time.time()
    # https://igraph.org/python/doc/api/igraph._igraph.GraphBase.html
    # community_walktrap through community_edge_betweenness
    # end_time = time.time()
    # print(end_time - start_time)
    # maintain cluster size, dynamic resolution?
    
    # https://leidenalg.readthedocs.io/en/stable/reference.html
    partition = la.find_partition(G, la.RBConfigurationVertexPartition, resolution_parameter=160)
    # leidenalg.find_partition(graph, partition_type, initial_membership=None, weights=None, n_iterations=2, max_comm_size=0, seed=None, **kwargs)

    membership = np.array(partition.membership)
    print(f'    Number of clusters: {len(np.unique(membership))}')
    print(f'    Modularity: {partition.modularity}')

    # Create the Leiden CSR graph
    print("    Leiden graph...")
    G.contract_vertices(membership, combine_attrs=None)
    G.simplify(multiple=False, loops=True, combine_edges=None)
    G.es["weight"] = [1 - w for w in G.es["weight"]] # cosine distance -> similarity
    G.simplify(multiple=True, loops=True, combine_edges=sum)
    i, j = zip(*G.get_edgelist())
    ic_edges = sp.csr_matrix((G.es["weight"], (i,j)), shape=(G.vcount(), G.vcount()))
    ic_edges = ic_edges + ic_edges.T
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
        return(knn_matrix2indist(C))
        
    knn_indices, knn_dists = my_knn(ic_edges)

    print("    Leiden UMAP...")
    from umap import UMAP        
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   random_state = None,
                   verbose = False,
                   low_memory=True,
                   n_epochs = 5000,
                   precomputed_knn = (knn_indices, knn_dists),
                   init = "spectral",
                   n_jobs = -1
                  )
    embedding = reducer.fit_transform(ic_edges)
    embedding -= np.mean(embedding, axis=0) # center

    return embedding, membership

embedding, membership = leiden_init(G)
init = embedding[membership]

print("Saving Leiden initialization...")

# Initialization plot
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(embedding[:, 0], embedding[:, 1], s=1)
ax.set_title('Leiden Initialization')
ax.axis('equal')
fig.savefig(os.path.join(out_dir, "leiden.pdf"), dpi=200)
del fig, ax

np.savez_compressed(os.path.join(out_dir, 'leiden.npz'), embedding=embedding, membership=membership)
# save G too? del G?

### Load + subset data #########################################################

# Load matrix
print('\nLoading the matrix...')
df = pd.read_csv(os.path.join(in_dir, 'matrix.csv.gz'), compression='gzip')
df.sb1_index -= 1 # convert from 1- to 0-indexed
df.sb2_index -= 1 # convert from 1- to 0-indexed
# assert set(df.sb1_index) == set(range(sb1.shape[0]))
# assert set(df.sb2_index) == set(range(sb2.shape[0]))

if bead == 1:
    mat = sp.csr_matrix((df['umi'], (df['sb1_index'], df['sb2_index'])))
if bead == 2:
    mat = sp.csr_matrix((df['umi'], (df['sb2_index'], df['sb1_index'])))
del df ; gc.collect()

# Subset matrix
assert len(knn_mask.mask) == mat.shape[0]
assert type(mat) == sp._csr.csr_matrix
mat = mat[knn_mask.mask]
print(f"Final matrix dimension: {mat.shape}")
print(f"Final matrix size: {csr_nbytes(mat)/1024/1024:.2f} MiB")
gc.collect()

# Fill in -1 KNN elements with neighbors of neighbors
# from sklearn.metrics.pairwise import cosine_similarity
# b_list = np.where(np.any(knn_indices == -1, axis=1))[0]
# print(f"Filling out neighbors for {len(b_list)} beads")
# for b in b_list:
#     ns = G.neighborhood(b, order=2, mindist=2)
#     cosine_similarity(mat[b], mat[ns])

# Load barcodes
print('Loading the barcodes...')
sb = pd.read_csv(os.path.join(in_dir, f"sb{bead}.txt.gz"), header=None, names=['sb'])

# Subset barcodes
assert len(knn_mask.mask) == sb.shape[0]
assert type(sb) == pd.core.frame.DataFrame
sb = sb[knn_mask.mask]
gc.collect()

assert mat.shape[0] == knn_indices.shape[0] == knn_dists.shape[0] == sb.shape[0]

### UMAP TIME ##################################################################

print("\nRunning UMAP...")
embeddings = []
if opts["n_epochs"] <= 1000:
    embeddings.append(my_umap(mat, (knn_indices, knn_dists), init, opts, cores))
else:
    n_epochs = opts["n_epochs"]
    opts["n_epochs"] = 1000
    embeddings.append(my_umap(mat, (knn_indices, knn_dists), init, opts, cores))
    for i in range(int(np.ceil(n_epochs/1000))-1):
        embeddings.append(my_umap(mat, (knn_indices, knn_dists), embeddings[-1], opts, cores))
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
names = ["umap", "convergence", "umaps", "umap_umi", "leiden", "knn"]
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
