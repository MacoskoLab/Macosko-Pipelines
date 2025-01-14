import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

### PLOTTING METHODS ###########################################################

# (x, y)
def hexmap(embedding, weights=None, title="", fontsize=12):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    x, y = embedding[:, 0], embedding[:, 1]
    hb = ax.hexbin(x, y, cmap='viridis' if weights is None else 'inferno', C = weights, linewidths=0.5)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{title}', fontsize=fontsize)
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.tight_layout()
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

# # TODO: highlight reciprocated ones
# def embedding_neighbor_locations(knn_indices, knn_dists, embedding, nn=45, n=16):
#     from matplotlib.cm import viridis
#     from matplotlib.colors import Normalize
#     assert knn_indices.shape[0] == knn_dists.shape[0] == embedding.shape[0]
#     assert knn_indices.shape[1] == knn_dists.shape[1]
    
#     # Create the grid
#     nrows = np.ceil(np.sqrt(n)).astype(int)
#     fig, axes = plt.subplots(nrows, nrows, figsize=(8, 8))
#     if type(axes) != np.ndarray:
#         axes = np.array([axes])
#     axes = axes.flatten()

#     selected_indices = np.random.randint(0, embedding.shape[0], size=n)
#     x = embedding[:,0] ; y = embedding[:,1]

#     # Plot the data
#     for ax, i in zip(axes, selected_indices):
#         indices = knn_indices[i, 1:nn]
#         dists = knn_dists[i, 1:nn]
#         colors = viridis(Normalize(vmin=0, vmax=1)(dists))

#         ax.scatter(x, y, color='grey', s=10, alpha=0.5)
#         ax.scatter(x[indices], y[indices], color=colors, s=20)
#         ax.scatter(x[i], y[i], color='red', s=30)

#         ax.set_xlim(min(x[indices]), max(x[indices]))
#         ax.set_ylim(min(y[indices]), max(y[indices]))
#         ax.set_aspect('equal', adjustable='box')
#         ax.set_title(i)
    
#     for ax in axes[n:]:
#         ax.set_visible(False)
    
#     fig.tight_layout()
#     return fig, axes

# def embedding_neighbor_distances(knn_indices, knn_dists, embedding, nn=45):
#     assert knn_indices.shape[0] == knn_dists.shape[0] == embedding.shape[0]
#     assert knn_indices.shape[1] == knn_dists.shape[1]

#     dists = np.vstack([np.linalg.norm(embedding[inds] - embedding[inds[0]], axis=1) for inds in knn_indices])
    
#     def n_hist(ax, dists, nn):
#         data = dists[:,nn]
#         ax.hist(np.log10(data), bins=100)
#         ax.set_xlabel('UMAP distance (log10)')
#         ax.set_ylabel('Count')
#         ax.set_title(f'UMAP Distance to neighbor {nn}')
#         meanval = np.log10(np.mean(data))
#         ax.axvline(meanval, color='red', linestyle='dashed')
#         ax.text(meanval+0.1, ax.get_ylim()[1] * 0.95, f'Mean: {10**meanval:.2f}', color='black', ha='left')

#     fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    
#     n_hist(axes[0,0], dists, 1)
#     n_hist(axes[0,1], dists, nn)
    
#     ax=axes[1,0]
#     data = np.mean(dists[:,1:nn], axis=1)
#     ax.hist(np.log10(data), bins=100)
#     ax.set_xlabel('UMAP distance (log10)')
#     ax.set_ylabel('Count')
#     ax.set_title(f'UMAP Distance to average neighbor (45)')
#     meanval = np.log10(np.mean(data))
#     ax.axvline(meanval, color='red', linestyle='dashed')
#     ax.text(meanval+0.1, ax.get_ylim()[1] * 0.95, f'Mean: {10**meanval:.2f}', color='black', ha='left')
    
#     axes[1,1].hexbin(np.log10(dists[:,1:nn]), knn_dists[:,1:nn], gridsize=100, bins='log', cmap='plasma')
#     axes[1,1].set_xlabel('UMAP distance (log10)')
#     axes[1,1].set_ylabel('Cosine Distance')
#     axes[1,1].set_title(f'Cosine vs. UMAP Distance ({nn})')
    
#     fig.tight_layout()
    
#     return fig, axes

# TODO: annotate number of beads with incomplete neighbors
def knn_plot(knn_indices, knn_dists, tlu):
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    def my_histogram(ax, vec, title, xlab):
        ax.hist(vec, bins=100)
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel("Counts")
        
        xmean = np.mean(vec)
        ax.axvline(xmean, color='red', linestyle='dotted', linewidth=1.5)
        ax.text(xmean + ax.get_xlim()[1] * 0.02, ax.get_ylim()[1] * 0.95,
                f'Mean: {xmean:.2f}', color='red', ha='left')

    # Nearest neighbor distance
    nnd = knn_dists[:,1]
    nnd = nnd[np.isfinite(nnd)]
    my_histogram(axs[0,0], nnd, "Nearest neighbor distance", "Cosine distance")
    # Furthest neighbor distance
    fnd = knn_dists[:,-1]
    fnd = fnd[np.isfinite(fnd)]
    my_histogram(axs[0,1], fnd, f"Furthest neighbor ({knn_dists.shape[1]}) distance", "Cosine distance")
    # Number of in-edges
    unique, counts = np.unique(knn_indices[:,1:], return_counts=True)
    my_histogram(axs[1,0], counts[unique >= 0], f"Number of in-edges", "In-edges")
    # Clustering coefficient
    my_histogram(axs[1,1], tlu, f"Clustering coefficient", "transitivity_local_undirected")

    plt.tight_layout()
    return fig, axs

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

### KNN METHODS ################################################################

class KNNMask:
    # initialize the valid indices
    def __init__(self, knn_matrix):
        assert knn_matrix.shape[0] == knn_matrix.shape[1]
        self.mask = np.ones(knn_matrix.shape[0], dtype=bool)

    # input knn_matrix and mask, return filtered matrix
    def apply_mask(self, knn_matrix, mask):
        assert knn_matrix.shape[0] == knn_matrix.shape[1] == len(mask)
        assert mask.dtype == bool
        
        # update mask
        v = np.where(self.mask)[0]
        assert len(v) == len(mask)
        self.mask[v[~mask]] = False
        
        # filter matrix
        return(knn_matrix[mask,:][:,mask])

# Convert knn_matrix to (knn_indices, knn_dists)
def knn_matrix2indist(knn_matrix):
    lil = knn_matrix.tolil(copy=False)
    nrows = lil.shape[0]
    ncols = max(len(row) for row in lil.rows) + 1
    
    knn_indices = np.full((nrows, ncols), -1, dtype=np.int32)
    knn_dists = np.full((nrows, ncols), np.inf, dtype=np.float32)
    
    for i in range(nrows):
        inds = np.array([i]+lil.rows[i], dtype=np.int32)
        vals = np.array([0]+lil.data[i], dtype=np.float32)

        sorted_indices = np.argsort(vals)
        inds = inds[sorted_indices]
        vals = vals[sorted_indices]

        knn_indices[i, :len(inds)] = inds
        knn_dists[i, :len(vals)] = vals

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

    return sp.csr_matrix((data, (rows, cols)), shape=csr.shape)

def csr_nbytes(csr):
    assert type(csr) == sp._csr.csr_matrix
    return(csr.data.nbytes + csr.indptr.nbytes + csr.indices.nbytes)

def validate_knn_indist(knn_indices, knn_dists):
    assert knn_indices.shape == knn_dists.shape
    assert knn_indices.dtype == np.int32 and knn_dists.dtype == np.float32
    assert np.array_equal(knn_indices[:,0], np.arange(len(knn_indices)))
    assert np.all(-1 <= knn_indices) and np.all(knn_indices < len(knn_indices))
    assert np.all(knn_dists[:,0] == 0)
    assert np.all(knn_dists[:,1:] > 0)
    assert not np.any(np.isnan(knn_dists))

# do checks

### UMAP METHODS ###############################################################

def my_umap(mat, knn, init, opts, n_jobs=-1):
    from umap import UMAP
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   random_state = None,
                   verbose = True,
                   low_memory=True,
                   
                   n_neighbors = opts["n_neighbors"],
                   min_dist = opts["min_dist"],
                   spread = opts["spread"],
                   local_connectivity = opts["local_connectivity"],
                   repulsion_strength = opts["repulsion_strength"],
                   negative_sample_rate = opts["negative_sample_rate"],
                   n_epochs = opts["n_epochs"],

                   precomputed_knn = knn,
                   init = init,
                   n_jobs = n_jobs
                  )
    embedding = reducer.fit_transform(np.log1p(mat))
    return(embedding)
