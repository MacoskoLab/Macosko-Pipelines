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
import umap.plot
from umap.umap_ import nearest_neighbors

def get_args():
    parser = argparse.ArgumentParser(description='process recon seq data')
    parser.add_argument("-i", "--in_dir", help="input data folder", type=str, default=".")
    parser.add_argument("-o", "--out_dir", help="output data folder", type=str, default=".")
    parser.add_argument("-gs", "--gspath", help="gcloud storage path to cache data output", type=str, default="")
    # parser.add_argument("-c", "--core", help="define core type to use (CPU or GPU)", type=str, default="CPU")
    # bead type? tags or seq? "-e", "--exptype", help="define experiment type (seq or tags)", type=str, required=True,
    parser.add_argument("-c1", "--cutoff1", help="R1 UMI cutoff", type=int, default=0)
    parser.add_argument("-c2", "--cutoff2", help="R2 UMI cutoff", type=int, default=0)

    parser.add_argument("-a", "--algorithm", help="dimensionality reduction algo", type=str, default="umap")
    
    parser.add_argument("-n", "--n_neighbors", help="the number of neighboring sample points used for manifold approximation", type=int, default=25)
    parser.add_argument("-d", "--min_dist", help="the effective minimum distance between embedded points", type=float, default=0.99)
    parser.add_argument("-s", "--spread", help="the effective scale of embedded points", type=float, default=1.0)
    parser.add_argument("-I", "--init", help="how to initialize the low dimensional embedding", type=str, default="spectral")
    parser.add_argument("-m", "--metric", help="the metric to use to compute distances in high dimensional space", type=str, default="cosine")
    parser.add_argument("-N", "--n_epochs", help="the number of training epochs to be used in optimizing the low dimensional embedding", type=int, default=10000)
    
    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

args = get_args()
in_dir = args.in_dir ; assert all(os.path.isfile(os.path.join(in_dir, file)) for file in ['matrix.csv.gz', 'sb1.txt.gz', 'sb2.txt.gz'])
c1 = args.cutoff1 ; print(f"R1 UMI cutoff = {c1}")
c2 = args.cutoff2 ; print(f"R2 UMI cutoff = {c2}")
base = f"ANCHOR_c1={c1}_c2={c2}"

algo = args.algorithm ; print(f"algorithm = {algo}")
if algo == "umap":
    n_neighbors = args.n_neighbors ; print(f"n_neighbors = {n_neighbors}")
    min_dist = args.min_dist       ; print(f"min_dist = {min_dist}")
    spread = args.spread           ; print(f"spread = {spread}")
    n_epochs = args.n_epochs       ; print(f"n_epochs = {n_epochs}")
    init = args.init               ; print(f"init = {init}")
    metric = args.metric           ; print(f"metric = {metric}")
    out_dir = os.path.join(args.out_dir, f"{base}_n={n_neighbors}_d={min_dist}_s={spread}_I={init}_m={metric}")
else:
    out_dir = os.path.join(args.out_dir, f"{base}")
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

# Filter the matrix
print("\nFiltering the beads...")
umi_before = sum(df["umi"])
if c1 > 0:
    grouped = df.groupby('sb1')['umi'].sum()
    sb1_keep = grouped[grouped <= c1].index
    print(f"{len(sb1)-len(sb1_keep)} R1 beads filtered ({round((len(sb1)-len(sb1_keep))/len(sb1)*100, 2)}%)")
if c2 > 0:
    grouped = df.groupby('sb2')['umi'].sum()
    sb2_keep = grouped[grouped <= c2].index
    print(f"{len(sb2)-len(sb2_keep)} R2 beads filtered ({round((len(sb2)-len(sb2_keep))/len(sb2)*100, 2)}%)")
if c1 > 0:
    df = df[df['sb1'].isin(sb1_keep)]
if c2 > 0:
    df = df[df['sb2'].isin(sb2_keep)]
if c1 > 0 or c2 > 0:
    codes1, uniques1 = pd.factorize(df['sb1'], sort=True)
    df['sb1'] = codes1
    codes2, uniques2 = pd.factorize(df['sb2'], sort=True)
    df['sb2'] = codes2
umi_after = sum(df["umi"])
print(f"{umi_before-umi_after} UMIs filtered ({round((umi_before-umi_after)/umi_before*100, 2)}%)")
assert sorted(list(set(df.sb1))) == list(range(len(set(df.sb1))))
assert sorted(list(set(df.sb2))) == list(range(len(set(df.sb2))))

# Rows are the anchor beads I wish to recon
# Columns are the features used for judging similarity
mat = coo_matrix((df['umi'], (df['sb2'], df['sb1'])))

# Get the previous embeddings
try:
    import gcsfs
    file_path = os.path.join(args.gspath, out_dir, "embeddings.npz")
    # file_path = r"gs://fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36/reconstruction/240615_SL-EXG_0144_A22KH5WLT3/D703_D704_D705_D706/ANCHOR_c1=0_c2=0_n=25_d=0.001_s=1.0_I=spectral_m=cosine/embeddings.npz"
    print(f"Searching {file_path}...")
    with gcsfs.GCSFileSystem().open(file_path, 'rb') as f:
        data = np.load(f)
        embeddings = [data[key] for key in data]
    print(f"{len(embeddings)} previous embeddings found")
except Exception as e:
    embeddings = []
    print(f"Embeddings load error: {str(e)}")
    print("No previous embeddings found, starting from scratch")

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
    if len(embeddings) == 0:
        embeddings.append(my_umap(mat, n_epochs=10))
        embeddings.append(my_umap(mat, n_epochs=90, init=embeddings[-1]))
        embeddings.append(my_umap(mat, n_epochs=900, init=embeddings[-1]))
    else:
        embeddings.append(my_umap(mat, n_epochs=1000))
    
    for i in range(ceil(n_epochs/1000)-1):
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
    sbs = [sb2[i] for i in uniques2] if (c1 > 0 or c2 > 0) else sb2
    assert embeddings[-1].shape[0] == len(sbs)
    with open(os.path.join(out_dir, "Puck.csv"), mode='w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(sbs)):
            writer.writerow([sbs[i], embeddings[-1][i,0], embeddings[-1][i,1]])
    
    print("\nDone!")
