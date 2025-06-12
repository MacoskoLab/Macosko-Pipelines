import os
import gc
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pypdf import PdfWriter
from helpers import *

def get_args():
    parser = argparse.ArgumentParser(description='2D embedding of reconstruction data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-c", "--cores", type=int, default=-1)
    parser.add_argument("-b", "--bead", type=int, default=2)
    parser.add_argument("-D", "--diameter", type=float, default=0) #in microns, <= 0 estimates using bead count
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
diameter = args.diameter         ; print(f"diameter: {'auto' if diameter <= 0 else diameter}")
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
name += (f"_D{diameter}" if diameter > 0 else "")
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

# Load sb
sb = pd.read_csv(os.path.join(in_dir, f"sb{bead}.txt.gz"), header=None, names=['sb'])['sb'].to_numpy()

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

# Compute scaling factor
if diameter <= 0:
    diameter_micron = estimate_diameter(in_dir)
    print(f"Estimated puck diameter: {diameter_micron}")
else:
    diameter_micron = diameter

# Create KNN object
knn = KNN(knn_matrix, sb, umi, opts["n_neighbors"])
del knn_matrix, sb, umi ; gc.collect()

# Initial bead filter
knn.filter_toofewneighbors()
knn.filter_disconnected()

### Load KNN graph ############################################################

print("\nCreating KNN graph...")
G = knn.graph()

print("\nComputing KNN metrics...")
knn_indices, knn_dists = knn.inddist()
nnd = knn_dists[:,0] ; nnd = nnd[np.isfinite(nnd)]
fnd = knn_dists[:,-1] ; fnd = fnd[np.isfinite(fnd)]
ies = np.array(G.degree(mode="in"))
tlu = np.array(G.transitivity_local_undirected(mode="zero", weights='weight'))
sts = np.array(G.strength(weights="weight"))
con = np.array(G.constraint(weights="weight"))
del knn_indices, knn_dists

print("\nPlotting KNN metrics...")
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
del fig, axs, my_histogram

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

# Intermediate bead filter
knn.filter_toofewneighbors()
knn.filter_disconnected()

### Leiden Initialization ######################################################

print("\nGenerating Leiden embedding...")
assert G.vcount() == knn.matrix.shape[0] == knn.matrix.shape[1]
membership, ic_edges, embedding = leiden_init(G)
assert G.vcount() == embedding.shape[0] == np.unique(membership).shape[0]
assert len(knn.membership) == len(membership)
knn.membership = membership
del G, membership ; gc.collect()

# Plot Leiden initialization
save_leiden_plots(embedding, knn, ic_edges, out_dir)

print("\nLeiden filter...")
# TODO

# Final bead filter
knn.filter_toofewneighbors()
knn.filter_disconnected()

# NO MORE FILTERING PAST THIS POINT!

# Plot filtering results
fig, ax = plt.subplots(figsize=(8, 8))
top_lines = [in_dir, name, f"Estimated diameter: {diameter_micron}µm"]
bot_lines = ["<Beads filtered>"] + [f"{n}: {b}" for n,b in knn.history]
lines = top_lines + [""] + bot_lines
[ax.text(0.01, 1-0.05*i, line, fontsize=10, va='top', ha='left') for i,line in enumerate(lines)]
ax.axis('off')
fig.savefig(os.path.join(out_dir, "name.pdf"))
del fig, ax, lines, top_lines, bot_lines

print("\nInitializing bead coordinates...")
# Generate bead initialization - stacked
'''
init = embedding[knn.membership]
'''
# Generate bead initialization - random jitter (default)
'''
init = embedding[knn.membership]
from scipy.spatial import KDTree
sd = np.mean(KDTree(embedding).query(embedding, k=2)[0][:,1]) / 1000
init += np.random.normal(loc=0, scale=sd, size=(init.shape[0],2))
'''
# Generate bead initialization - weighted average
knn_indices, knn_dists = knn.inddist()
weights = 1-knn_dists ; weights[knn_indices < 0] = 0
init = embedding[knn.membership]
for _ in range(3):
    init = np.stack((np.average(init[knn_indices, 0], axis=1, weights=weights),
                     np.average(init[knn_indices, 1], axis=1, weights=weights)), axis=1)
del knn_indices, knn_dists, weights

print("\nSaving initialization...")
np.savez_compressed(os.path.join(out_dir, 'init.npz'), embedding=embedding, membership=knn.membership, init=init, ic_edges=ic_edges)
assert init.shape[0] == knn.matrix.shape[0] == knn.matrix.shape[1]
del embedding, ic_edges ; gc.collect()

### UMAP TIME ##################################################################

print("\nRunning UMAP...")
from umap import UMAP
reducer = UMAP(n_components = 2,
               metric = "precomputed",
               random_state = None,
               verbose = True,
               init = init,
               n_jobs = cores,
               **opts
              )
embedding = reducer.fit_transform(knn.matrix.maximum(knn.matrix.T))
assert np.all(np.isfinite(embedding))
embedding -= np.mean(embedding, axis=0)

### WRITE RESULTS ##############################################################

print("\nWriting results...")
knn.check()
assert embedding.shape[0] == knn.matrix.shape[0]

# Create the Puck file
embedding_scaled = embedding / np.mean(np.ptp(embedding, axis=0)) * diameter_micron
np.savetxt(os.path.join(out_dir, 'Puck.csv'),
           np.column_stack([knn.sb, embedding_scaled]),
           delimiter=',', fmt='%s')

# Plot the UMAP embedding
title = f"UMAP hexbin ({embedding.shape[0]:} beads)"
fig, ax = hexmap(embedding, title=title)
fig.savefig(os.path.join(out_dir, "umap.pdf"), dpi=200)

# Plot the UMAP embedding (UMI-weighted)
title = f"UMAP hexbin (UMI-weighted)"
fig, ax = hexmap(embedding, knn.umi, title=title, legend=True)
fig.savefig(os.path.join(out_dir, "umap_umi.pdf"), dpi=200)

# Two additional informational plots
save_umap_neighbor_plots(embedding_scaled, knn, out_dir)

# Make summary pdf
# creation order: knn, init, init2, name, umap, umap_umi, umap_stats, umap_neighborhoods
names = ["umap", "umap_umi", "knn", "init", "init2", "umap_dists", "umap_neighbors", "name"]
paths = [os.path.join(out_dir, n+".pdf") for n in names]
files = [p for p in paths if os.path.isfile(p)]
if len(files) > 0:
    merger = PdfWriter()
    [merger.append(file_name) for file_name in files]
    merger.write(os.path.join(out_dir, "summary.pdf"))
    merger.close()
    [os.remove(file) for file in files]

print("Done!")
