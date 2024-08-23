import os
import gc
import sys
import csv
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from PyPDF2 import PdfMerger
from helpers import *

#os.chdir("/home/nsachdev/recon/data/240609_SL-EXC_0308_A22KHFYLT3__240716_SL-EXG_0167_B22CG3TLT4/D702-1234_D703_D706_D707-5__D702-8_D703-2_D706-2_D707-8/")
# os.chdir("/home/nsachdev/recon/data/240609_SL-EXC_0308_A22KHFYLT3/D707-1234")
os.chdir("/home/nsachdev/recon/data/240411/H14_1")
# os.chdir("/home/nsachdev/recon/data/240615_SL-EXG_0144_A22KH5WLT3/D703_D704_D705_D706")
# os.chdir("/home/nsachdev/recon/data/240609_SL-EXC_0308_A22KHFYLT3/D702-1234_D703_D706_D707-5")

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-gs", "--gspath", help="gcloud storage path to cache data output", type=str, default="")
    parser.add_argument("-u", "--unit", help="define core type to use (CPU or GPU)", type=str, default="CPU")

    parser.add_argument("-n", "--n_neighbors", help="the number of neighboring points used for manifold approximation", type=int, default=45)
    parser.add_argument("-d", "--min_dist", help="the effective minimum distance between embedded points", type=float, default=0.1)
    parser.add_argument("-I", "--init", help="how to initialize the low dimensional embedding", type=str, default="spectral")
    parser.add_argument("-N", "--n_epochs", help="the number of epochs to be used in optimizing the embedding", type=int, default=5000)
    parser.add_argument("-c", "--connectivity", help="'none', 'min_tree', or 'full_tree'", type=str, default="full_tree")
    parser.add_argument("-n2", "--n_neighbors2", help="the new NN to pick for MNN", type=int, default=45)
    
    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

args = get_args()

in_dir = args.in_dir             ; print(f"input directory = {in_dir}")
out_dir = args.out_dir           ; print(f"output directory = {out_dir}")
gspath = args.gspath             ; print(f"gs:// output path = {gspath}")
unit = args.unit                 ; print(f"processing unit = {unit}")

n_neighbors = args.n_neighbors   ; print(f"n_neighbors = {n_neighbors}")
min_dist = args.min_dist         ; print(f"min_dist = {min_dist}")
init = args.init                 ; print(f"init = {init}")
n_epochs = args.n_epochs         ; print(f"n_epochs = {n_epochs}")
connectivity = args.connectivity ; print(f"connectivity = {connectivity}")
n_neighbors2 = args.n_neighbors2 ; print(f"n_neighbors2 = {n_neighbors2}")

name = f"UMAP_n={n_neighbors}_d={min_dist}_I={init}"

if connectivity != "none":
    assert connectivity in ["min_tree", "full_tree"]
    name += f"_c={connectivity.replace('_', '')}{n_neighbors2}"
    assert n_neighbors2 <= n_neighbors

print(f"name = {name}")
out_dir = os.path.join(out_dir, name)

assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.csv.gz', 'sb2.csv.gz'])
os.makedirs(out_dir, exist_ok=True)
assert os.path.exists(out_dir)
print(f"output directory = {out_dir}")

# Get the previous embeddings
print("\nDownloading previous embeddings...")
file_path = os.path.join(gspath, name, "embeddings.npz")
print(f"Searching {file_path}...")
try:
    import gcsfs
    with gcsfs.GCSFileSystem().open(file_path, 'rb') as f:
        data = np.load(f)
        embeddings = [data[key] for key in data]
    print(f"{len(embeddings)} previous embeddings found")
except Exception as e:
    embeddings = []
    print(f"Embeddings load error: {str(e)}")
    print("No previous embeddings found, starting from scratch")

metadata = {}

sys.stdout.flush()

### Load the data ##############################################################

print("\nReading the matrix...")
df = pd.read_csv(os.path.join(in_dir, 'matrix.csv.gz'), compression='gzip')
df.sb1_index -= 1 # convert from 1- to 0-indexed
df.sb2_index -= 1 # convert from 1- to 0-indexed
sb1 = pd.read_csv(os.path.join(in_dir, 'sb1.csv.gz'), compression='gzip')
sb2 = pd.read_csv(os.path.join(in_dir, 'sb2.csv.gz'), compression='gzip')
assert sorted(list(set(df.sb1_index))) == list(range(sb1.shape[0]))
assert sorted(list(set(df.sb2_index))) == list(range(sb2.shape[0]))

metadata["init"] = {"sb1": sb1.shape[0], "sb2": sb2.shape[0], "umi": sum(df["umi"])}
fig, axes = uvc(df) ; fig.savefig(os.path.join(out_dir,'uvc.pdf'), format='pdf') ; del fig, axes
print(f"{sb1.shape[0]} R1 barcodes")
print(f"{sb2.shape[0]} R2 barcodes")

print("\nFiltering the beads...")
df, uniques1, uniques2, fig, meta = connection_filter(df)
fig.savefig(os.path.join(out_dir,'connections.pdf'), format='pdf') ; del fig
metadata["connection_filter"] = meta ; del meta

# Rows are the beads you wish to recon
# Columns are the features used for judging similarity
mat = coo_matrix((df['umi'], (df['sb2_index'], df['sb1_index']))).tocsr()
del df

sys.stdout.flush()

### Compute the KNN ############################################################

print("\nComputing the KNN...")
knn_indices1, knn_dists1 = knn_descent(np.log1p(mat), n_neighbors)

print("\nFiltering the KNN...")
filter_indexes, fig, meta = knn_filter(knn_indices1, knn_dists1)
fig.savefig(os.path.join(out_dir,'knn.pdf'), format='pdf') ; del fig
metadata["knn_filter"] = meta ; del meta
m = np.array([i not in filter_indexes for i in np.arange(mat.shape[0])], dtype=bool)
mat = mat[m,:] ; uniques2 = uniques2[m]

print("\nRe-computing the KNN...")
knn_indices2, knn_dists2 = knn_descent(np.log1p(mat), n_neighbors)
knn_indices, knn_dists = knn_merge(knn_indices1[m,:], knn_dists1[m,:], knn_indices2, knn_dists2)
np.savez_compressed(os.path.join(out_dir, "knn.npz"), indices=knn_indices, dists=knn_dists)
knn = (knn_indices, knn_dists)

if connectivity != "none":
    mnn_indices, mnn_dists = mutual_nn_nearest(knn_indices, knn_dists, n_neighbors, n_neighbors2, connectivity)
    np.savez_compressed(os.path.join(out_dir, "mnn.npz"), indices=mnn_indices, dists=mnn_dists)
    assert np.all(np.isfinite(mnn_indices))
    assert np.all(np.isfinite(mnn_dists))
    n_neighbors = n_neighbors2
    knn = (mnn_indices, mnn_dists)

print(f"Final matrix dimension: {mat.shape}")
print(f"Final matrix size: {mat.data.nbytes/1024/1024:.2f} MiB")

### UMAP TIME ##################################################################

if unit.upper() == "CPU":
    from umap import UMAP
elif unit.upper() == "GPU":
    from cuml.manifold.umap import UMAP
else:
    exit(f"Unrecognized --unit flag {unit}")



def cpu_umap(mat, knn, n_epochs, init=init):
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   spread = 1.0,
                   random_state = None,
                   verbose = True,
                   precomputed_knn = knn,
                   
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   n_epochs = n_epochs,
                   init = init
                  )
    embedding = reducer.fit_transform(np.log1p(mat))
    return(embedding)

print("\nRunning UMAP...")
embeddings.append(my_umap(mat, knn, n_epochs=1000))    
for i in range(int(np.ceil(n_epochs/1000))-1):
    print(i+2)
    # # Upload intermediate embeddings
    # try:
    #     import gcsfs
    #     file_path = os.path.join(gspath, name, "embeddings.npz")
    #     with gcsfs.GCSFileSystem().open(file_path, 'wb') as f:
    #         np.savez(f, **{f"arr_{i}":e for i,e in enumerate(embeddings)})
    #     print("Intermediate embeddings successfully uploaded")
    # except Exception as e:
    #     print(f"Unable to upload intermediate embeddings: {str(e)}")
    embeddings.append(my_umap(mat, knn, init=embeddings[-1], n_epochs=1000))

print("\nWriting results...")
embedding = embeddings[-1]

# Save the embeddings
np.savez_compressed(os.path.join(out_dir, "embeddings.npz"), *embeddings)

# Create the Puck file
sbs = [sb2["sb2"][i] for i in uniques2]
assert embedding.shape[0] == len(sbs)
with open(os.path.join(out_dir, "Puck.csv"), mode='w', newline='') as file:
    writer = csv.writer(file)
    for i in range(len(sbs)):
        writer.writerow([sbs[i], embedding[i,0], embedding[i,1]])

# Plot the umap
title = f"umap hexbin ({embedding.shape[0]:} anchor beads) [{(len(embeddings))*1000} epochs]"
fig, ax = hexmap(embedding, title)
fig.savefig(os.path.join(out_dir, "umap.pdf"), dpi=200) ; del fig

# Plot the intermediate embeddings
fig, axes = hexmaps(embeddings, titles=[(i+1)*1000 for i in range(len(embeddings))])
fig.savefig(os.path.join(out_dir, "umaps.pdf"), dpi=200) ; del fig

# Plot the convergence
if len(embeddings) > 1:
    fig.savefig(os.path.join(out_dir, "convergence.pdf"), dpi=200) ; del fig
embeddings

# Make summary pdf
names = ["umap", "connections", "knn", "uvc", "convergence", "umaps"]
paths = [os.path.join(out_dir, n+".pdf") for n in names]
files = [p for p in paths if os.path.isfile(p)]
if len(files) > 0:
    merger = PdfMerger()

    for file_name in files:
        merger.append(file_name)

    merger.write(os.path.join(out_dir, "summary.pdf"))
    merger.close()
    [os.remove(file) for file in files]

print("\nDone!")
