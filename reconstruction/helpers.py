import os
import gc
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

### PLOTTING METHODS ###########################################################

# (x, y)
def hexmap(embedding, weights=None, title="", fontsize=12, legend=False):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    x, y = embedding[:, 0], embedding[:, 1]
    hb = ax.hexbin(x, y, cmap='viridis' if weights is None else 'inferno', C = weights, linewidths=0.5)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{title}', fontsize=fontsize)
    [spine.set_visible(False) for spine in ax.spines.values()]
    plt.tight_layout()

    if legend:
        cbar = fig.colorbar(hb, ax=ax, shrink=0.75)
    
    return fig, ax

# [(x, y)]
def hexmaps(embeddings, titles=[], fontsize=10):
    assert type(embeddings) == type(titles) == list
    n = np.ceil(len(embeddings)**0.5).astype(int)
    fig, axes = plt.subplots(n, n, figsize=(8, 8))
    if type(axes) != np.ndarray:
        axes = np.array([axes])
    axes = axes.flatten()

    if len(titles) == 0:
        titles = ["" for _ in range(len(embeddings))]
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)
    assert len(titles) == len(embeddings)
    
    for ax, embedding, title in zip(axes, embeddings, titles):
        x, y = embedding[:, 0], embedding[:, 1]
        hb = ax.hexbin(x, y, cmap='viridis', linewidths=0.5)
        ax.axis('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'{title}', fontsize=fontsize)
        for spine in ax.spines.values():
            spine.set_visible(False)
        
    for ax in axes[len(embeddings):]:
        ax.set_visible(False)

    return fig, axes

# (x, y, color)
def beadplot(embedding, colors, cmap='viridis'):
    assert embedding.shape[0] == len(colors)
    if isinstance(embedding, np.ndarray):
        puck = np.hstack((embedding, colors.reshape(-1,1)))
        puck = puck[puck[:,2].argsort()]
        x = puck[:,0]
        y = puck[:,1]
        c = puck[:,2]
    elif isinstance(embedding, pd.DataFrame):
        puck = embedding
        puck.loc[:,"color"] = colors
        puck = puck.sort_values(by=puck.columns[2])
        x = puck.iloc[:,0]
        y = puck.iloc[:,1]
        c = puck.iloc[:,2]
    else:
        raise TypeError("Input must be a NumPy ndarray or a pandas DataFrame")

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, c=c, cmap=cmap, s=0.1)
    plt.colorbar()
    plt.xlabel('xcoord')
    plt.ylabel('ycoord')
    plt.axis('square')
    return plt

# [embedding]
def convergence_plot(embeddings):
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    x = [(i+2)*1000 for i in range(len(embeddings)-1)]

    assert type(embeddings) == list
    if len(embeddings) < 2:
        return fig, axes

    y1 = [L2_distance(embeddings[i], embeddings[i+1]) for i in range(len(embeddings)-1)]
    axes[0].scatter(x, y1, color='blue')
    axes[0].plot(x, y1, color='red')
    axes[0].set_xlabel('Epochs')
    axes[0].set_ylabel('Mean UMAP Distance')
    axes[0].set_title('Average L2 Displacement')

    y2 = [procrustes_distance(embeddings[i], embeddings[i+1]) for i in range(len(embeddings)-1)]
    axes[1].scatter(x, y2, color='blue')
    axes[1].plot(x, y2, color='red')
    axes[1].set_xlabel('Epochs')
    axes[1].set_ylabel('Disparity')
    axes[1].set_title('Procrustes Disparity')
    
    plt.tight_layout()
    return fig, axes

### MISC. HELPERS ##############################################################

# embedding1, embedding2
def L2_distance(p1, p2):
    assert p1.shape == p2.shape
    dists = np.sqrt(np.sum(np.square(p1 - p2), axis=1))
    return np.sum(dists) / p1.shape[0]

# embedding1, embedding2
def procrustes_distance(p1, p2):
    assert p1.shape == p2.shape
    from scipy.spatial import procrustes
    _, _, disparity = procrustes(p1, p2)
    return disparity

def estimate_diameter(in_dir):
    with gzip.open(os.path.join(in_dir, 'sb1.txt.gz'), 'rt') as f:
        num_sb1 = sum(1 for _ in f)
    with gzip.open(os.path.join(in_dir, 'sb2.txt.gz'), 'rt') as f:
        num_sb2 = sum(1 for _ in f)
    n = num_sb1+num_sb2
    
    if n > 17*1e6:
        diameter_micron = 70_000
    elif n > 8.43*1e6:
        diameter_micron = 40_000
    elif n > 4.22*1e6:
        diameter_micron = 30_000
    elif n > 1.63*1e6:
        diameter_micron = 20_000
    elif n > 0.8*1e6:
        diameter_micron = 12_000
    else:
        sys.exit("Unknown puck size")

    return diameter_micron

### KNN METHODS ################################################################

class KNN:
    def __init__(self, knn_matrix, sb, umi):
        assert knn_matrix.shape[0] == knn_matrix.shape[1] == sb.shape[0] == umi.shape[0]
        self.matrix = knn_matrix
        self.sb = sb
        self.umi = umi
        self.mask = np.ones(knn_matrix.shape[0], dtype=bool)
        self.history = []

    # subset data to mask
    def apply_mask(self, mask, name):
        assert self.matrix.shape[0] == self.matrix.shape[1] == self.sb.shape[0] == self.umi.shape[0] == mask.shape[0]
        assert mask.dtype == bool and mask.ndim == 1

        # update data
        self.matrix = self.matrix[mask,:][:,mask]
        self.sb = self.sb[mask]
        self.umi = self.umi[mask]
        assert self.matrix.shape[0] == self.matrix.shape[1] == self.sb.shape[0] == self.umi.shape[0] == np.sum(mask)
        
        # update mask
        v = np.where(self.mask)[0]
        assert len(v) == len(mask)
        self.mask[v[~mask]] = False
        assert np.sum(self.mask) == np.sum(mask)

        # update history
        self.history.append((name, np.sum(~mask)))

        print(f"{np.sum(~mask)} beads removed")
        gc.collect()
        return None

# Convert knn_matrix to (knn_indices, knn_dists)
def knn_matrix2indist(knn_matrix, k=None):
    lil = knn_matrix.tolil(copy=False)
    nrows = lil.shape[0]
    ncols = max(len(row) for row in lil.rows) + 1 if (k is None) else k
    
    knn_indices = np.full((nrows, ncols), -1, dtype=np.int32)
    knn_dists = np.full((nrows, ncols), np.inf, dtype=np.float32)
    
    for i in range(nrows):
        inds = np.array([i]+lil.rows[i], dtype=np.int32)
        vals = np.array([0]+lil.data[i], dtype=np.float32)

        sorted_indices = np.argsort(vals)
        inds = inds[sorted_indices][:ncols]
        vals = vals[sorted_indices][:ncols]

        knn_indices[i, :len(inds)] = inds
        knn_dists[i, :len(vals)] = vals

    assert knn_indices.shape[0] == knn_matrix.shape[0]
    validate_knn_indist(knn_indices, knn_dists)
    return knn_indices, knn_dists

# convert (knn_indices, knn_dists) back into knn_matrix
def knn_indist2matrix(knn_indices, knn_dists):
    validate_knn_indist(knn_indices, knn_dists)
    rows = np.repeat(knn_indices[:,0], knn_indices.shape[1]-1)
    cols = knn_indices[:,1:].ravel()
    vals = knn_dists[:,1:].ravel()
    
    # remove missing elements before constructing matrix
    remove = (cols < 0) | (vals <= 0) | ~np.isfinite(cols) | ~np.isfinite(vals)
    rows = rows[~remove]
    cols = cols[~remove]
    vals = vals[~remove]
    if np.sum(remove) > 0:
        print(f"{np.sum(remove)} values removed during matrix construction")
    
    knn_matrix = sp.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
    
    return knn_matrix

# subset a knn_matrix to the closest k neighbors in each row
def csr_k_nearest(csr, k):
    assert type(csr) == sp._csr.csr_matrix
    rows, cols, data = [], [], []
    
    for i in range(csr.shape[0]):
        row_start, row_end = csr.indptr[i], csr.indptr[i+1]
        row_data = csr.data[row_start:row_end]
        row_indices = csr.indices[row_start:row_end]

        if len(row_data) > k:
            top_k_indices = np.argpartition(row_data, k)[:k]
            row_data = row_data[top_k_indices]
            row_indices = row_indices[top_k_indices]

        rows.extend([i] * len(row_data))
        cols.extend(row_indices)
        data.extend(row_data)

    return sp.csr_matrix((data, (rows, cols)), shape=csr.shape) # , dtype=np.float32

def validate_knn_indist(knn_indices, knn_dists):
    assert knn_indices.shape == knn_dists.shape
    assert knn_indices.dtype == np.int32 and knn_dists.dtype == np.float32
    assert np.array_equal(knn_indices[:,0], np.arange(len(knn_indices)))
    assert np.all(-1 <= knn_indices) and np.all(knn_indices < len(knn_indices))
    assert np.all(knn_dists[:,0] == 0)
    assert np.all(knn_dists[:,1:] > 0)
    assert not np.any(np.isnan(knn_dists))

# do checks
