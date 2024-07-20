import os
import gc
import sys
import csv
import gzip
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil
from scipy.sparse import coo_matrix
from umap import UMAP
from helpers import *

# os.chdir("/home/nsachdev/recon/data/6mm")
# os.chdir("/home/nsachdev/recon/data/1.2cm")
# os.chdir("/home/nsachdev/recon/data/2cm")
# os.chdir("/home/nsachdev/recon/data/3cm")

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-gs", "--gspath", help="gcloud storage path to cache data output", type=str, default="")
    # parser.add_argument("-c", "--core", help="define core type to use (CPU or GPU)", type=str, default="CPU")
    # tags or seq? "-e", "--exptype", help="define experiment type (seq or tags)", type=str, required=True,
    
    parser.add_argument("-l1", "--low1", help="R1 connection minimum", type=int, default=10)
    parser.add_argument("-l2", "--low2", help="R2 connection minimum", type=int, default=10)
    parser.add_argument("-h1", "--high1", help="R1 connection maximum", type=int, default=1000)
    parser.add_argument("-h2", "--high2", help="R2 connection maximum", type=int, default=1000)
    
    parser.add_argument("-a", "--algorithm", help="dimensionality reduction algo", type=str, default="UMAP") # UMAP
    
    parser.add_argument("-n", "--n_neighbors", help="the number of neighboring points used for manifold approximation", type=int, default=45)
    parser.add_argument("-d", "--min_dist", help="the effective minimum distance between embedded points", type=float, default=0.2)
    parser.add_argument("-s", "--spread", help="the effective scale of embedded points", type=float, default=1.0)
    parser.add_argument("-I", "--init", help="how to initialize the low dimensional embedding", type=str, default="spectral")
    parser.add_argument("-N", "--n_epochs", help="the number of epochs to be used in optimizing the embedding", type=int, default=5000)
    parser.add_argument("-c", "--connectivity", help="'none', 'min_tree', or 'full_tree'", type=str, default="none")
    parser.add_argument("-n2", "--n_neighbors2", help="the new NN to pick for MNN", type=int, default=45)
    
    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

args = get_args()
algo = args.algorithm ; print(f"algorithm = {algo}")
name = f"{algo}"

if algo == "UMAP":
    n_neighbors = args.n_neighbors ; print(f"n_neighbors = {n_neighbors}")
    min_dist = args.min_dist       ; print(f"min_dist = {min_dist}")
    spread = args.spread           ; print(f"spread = {spread}")
    init = args.init               ; print(f"init = {init}")
    n_epochs = args.n_epochs       ; print(f"n_epochs = {n_epochs}")
    name += f"_n={n_neighbors}_d={min_dist}_s={spread}_I={init}"
    
    connectivity = args.connectivity ; print(f"connectivity = {connectivity}")
    assert connectivity in ["none", "nearest", "min_tree", "full_tree"]
    if connectivity != "none":
        n_neighbors2 = args.n_neighbors2 ; print(f"n_neighbors2 = {n_neighbors2}")
        name += f"_c={connectivity.replace('_', '')}{n_neighbors2}"

l1 = args.low1  ; print(f"R1 connection minimum = {l1}")
l2 = args.low2  ; print(f"R2 connection minimum = {l2}")
h1 = args.high1 ; print(f"R1 connection maximum = {h1}")
h2 = args.high2 ; print(f"R2 connection maximum = {h2}")
name += f"_c1={l1}-{h1}_c2={l2}-{h2}"

print(f"name = {name}")

in_dir = args.in_dir
assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.csv.gz', 'sb2.csv.gz'])
print(f"input directory = {in_dir}")

out_dir = os.path.join(args.out_dir, name)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
assert os.path.exists(out_dir)
print(f"output directory = {out_dir}")

### Load the data ##############################################################

print("\nReading the matrix...")
df = pd.read_csv(os.path.join(in_dir, 'matrix.csv.gz'), compression='gzip', header=None, names=['sb1', 'sb2', 'umi'])
df.sb1 -= 1 # convert from 1- to 0-indexed
df.sb2 -= 1 # convert from 1- to 0-indexed
sb1 = pd.read_csv(os.path.join(in_dir, 'sb1.csv.gz'), compression='gzip', header=None, names=['sb1', 'umi', 'connections'])
sb2 = pd.read_csv(os.path.join(in_dir, 'sb2.csv.gz'), compression='gzip', header=None, names=['sb2', 'umi', 'connections'])
assert sorted(list(set(df.sb1))) == list(range(sb1.shape[0]))
assert sorted(list(set(df.sb2))) == list(range(sb2.shape[0]))
print(f"{sb1.shape[0]} R1 barcodes")
print(f"{sb2.shape[0]} R2 barcodes")

# Filter the matrix
print("\nFiltering the beads...")
umi_before = sum(df["umi"])
sb1_low  = np.where(sb1['connections'] <  l1)[0]
sb2_low  = np.where(sb2['connections'] <  l2)[0]
sb1_high = np.where(sb1['connections'] >= h1)[0]
sb2_high = np.where(sb2['connections'] >= h2)[0]
print(f"{len(sb1_low)} low R1 beads filtered ({len(sb1_low)/len(sb1)*100:.2f}%)")
print(f"{len(sb2_low)} low R2 beads filtered ({len(sb2_low)/len(sb2)*100:.2f}%)")
print(f"{len(sb1_high)} high R1 beads filtered ({len(sb1_high)/len(sb1)*100:.2f}%)")
print(f"{len(sb2_high)} high R2 beads filtered ({len(sb2_high)/len(sb2)*100:.2f}%)")
df = df[~df['sb1'].isin(sb1_low) & ~df['sb1'].isin(sb1_high) & ~df['sb2'].isin(sb2_low) & ~df['sb2'].isin(sb2_high)]
umi_after = sum(df["umi"])
print(f"{umi_before-umi_after} UMIs filtered ({(umi_before-umi_after)/umi_before*100:.2f}%)")
codes1, uniques1 = pd.factorize(df['sb1'], sort=True)
df.loc[:, 'sb1'] = codes1
codes2, uniques2 = pd.factorize(df['sb2'], sort=True)
df.loc[:, 'sb2'] = codes2
assert sorted(list(set(df.sb1))) == list(range(len(set(df.sb1))))
assert sorted(list(set(df.sb2))) == list(range(len(set(df.sb2))))

# Rows are the beads you wish to recon
# Columns are the features used for judging similarity
mat = coo_matrix((df['umi'], (df['sb2'], df['sb1']))).tocsr()
del df
print(f"Final matrix size: {mat.data.nbytes/1024/1024:.2f} MiB")
print(f"Final matrix dimension: {mat.shape}")

# # Get the previous embeddings
# print("\nDownloading previous embeddings...")
# file_path = os.path.join(args.gspath, name, "embeddings.npz")
# print(f"Searching {file_path}...")
# try:
#     import gcsfs
#     with gcsfs.GCSFileSystem().open(file_path, 'rb') as f:
#         data = np.load(f)
#         embeddings = [data[key] for key in data]
#     print(f"{len(embeddings)} previous embeddings found")
# except Exception as e:
#     embeddings = []
#     print(f"Embeddings load error: {str(e)}")
#     print("No previous embeddings found, starting from scratch")

embeddings = []
sys.stdout.flush()

### UMAP TIME ##################################################################

def my_umap(mat, n_epochs, init=init):
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   random_state = None,
                   low_memory = True,
                   verbose = True,
                   precomputed_knn = (knn_indices, knn_dists),
                   
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   spread = spread,
                   n_epochs = n_epochs,
                   init = init
                  )
    embedding = reducer.fit_transform(np.log1p(mat))
    return(embedding)

if algo == "UMAP":
    
    print("\nComputing the KNN...")
    knn_indices, knn_dists = knn_descent(np.log1p(mat), n_neighbors)
    
    if connectivity != "none":
        knn_indices, knn_dists = mutual_nn_nearest(knn_indices, knn_dists, n_neighbors, n_neighbors2, connectivity)
        assert np.all(np.isfinite(knn_indices))
        assert np.all(np.isfinite(knn_dists))
        n_neighbors = n_neighbors2
    
    print("\nRunning UMAP...")
    if len(embeddings) == 0:
        embeddings.append(my_umap(mat, n_epochs=20))
        embeddings.append(my_umap(mat, n_epochs=80, init=embeddings[-1]))
        embeddings.append(my_umap(mat, n_epochs=900, init=embeddings[-1]))
    else:
        embeddings.append(my_umap(mat, n_epochs=1000))
    
    for i in range(ceil(n_epochs/1000)-1):
        # # Upload intermediate embeddings
        # try:
        #     import gcsfs
        #     file_path = os.path.join(args.gspath, name, "embeddings.npz")
        #     with gcsfs.GCSFileSystem().open(file_path, 'wb') as f:
        #         np.savez(f, **{f"arr_{i}":e for i,e in enumerate(embeddings)})
        #     print("Intermediate embeddings successfully uploaded")
        # except Exception as e:
        #     print(f"Unable to upload intermediate embeddings: {str(e)}")

        # Run more umap
        print(i+2)
        embeddings.append(my_umap(mat, init=embeddings[-1], n_epochs=1000))

print("\nWriting results...")

# Save the embeddings
np.savez(os.path.join(out_dir, "embeddings.npz"), *embeddings)

# Plot the final UMAP
fig, ax = plt.subplots(figsize=(10, 8))
x, y = embeddings[-1][:, 0], embeddings[-1][:, 1]
hb = ax.hexbin(x, y, cmap='viridis', linewidths=0.1)
cb = fig.colorbar(hb, ax=ax, shrink = 0.75)
ax.set_title(f'umap hexbin ({embeddings[-1].shape[0]:,} anchor beads) [{(len(embeddings)-2)*1000} epochs]')
ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.axis('equal')
plt.tight_layout()
fig.savefig(os.path.join(out_dir, "umap.png"), dpi=200)
plt.close(fig)

# Create the Puck file
sbs = [sb2["sb2"][i] for i in uniques2]
assert embeddings[-1].shape[0] == len(sbs)
with open(os.path.join(out_dir, "Puck.csv"), mode='w', newline='') as file:
    writer = csv.writer(file)
    for i in range(len(sbs)):
        writer.writerow([sbs[i], embeddings[-1][i,0], embeddings[-1][i,1]])

print("\nDone!")
