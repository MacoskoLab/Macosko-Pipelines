import os
import gc
import sys
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
    elif n > 0.4*1e6:
        diameter_micron = 12_000
    else:
        sys.exit(f"Unknown puck size: sb1={num_sb1}, sb2={num_sb2}")

    return diameter_micron

### KNN METHODS ################################################################

class KNN:
    def __init__(self, matrix, sb, umi, n_neighbors):
        self.matrix = matrix
        self.sb = sb
        self.umi = umi
        self.membership = np.full(matrix.shape[0], -1)
        
        self.mask = np.ones(matrix.shape[0], dtype=bool)
        self.n_neighbors = n_neighbors
        self.history = []

        self.check()

    # Validate internal state
    def check(self):
        assert self.matrix.shape[0] == self.matrix.shape[1] == self.sb.shape[0]
        assert self.sb.shape[0] == self.umi.shape[0] == self.membership.shape[0]
        assert type(self.matrix) == sp._csr.csr_matrix and self.matrix.dtype == np.float32
        assert type(self.sb) == np.ndarray and self.sb.ndim == 1
        assert type(self.umi) == np.ndarray and self.umi.ndim == 1
        assert type(self.membership) == np.ndarray and self.membership.ndim == 1
        assert type(self.mask) == np.ndarray and self.mask.ndim == 1 and self.mask.dtype == bool
        assert np.diff(self.matrix.indptr).max() >= self.n_neighbors
        assert self.matrix.shape[0] > 0
    
    # Subset data to input mask
    def apply_mask(self, mask, name):
        assert type(mask) == np.ndarray and mask.dtype == bool and mask.ndim == 1
        assert self.matrix.shape[0] == mask.shape[0]
        self.check()

        # update data
        self.matrix = self.matrix[mask,:][:,mask]
        self.sb = self.sb[mask]
        self.umi = self.umi[mask]
        self.membership = self.membership[mask]
        assert self.matrix.shape[0] == np.sum(mask)
        
        # update mask
        v = np.where(self.mask)[0]
        assert len(v) == len(mask)
        self.mask[v[~mask]] = False
        assert np.sum(self.mask) == np.sum(mask)

        # update history
        self.history.append((name, np.sum(~mask)))

        print(f"{np.sum(~mask)} beads removed")
        self.check() ; gc.collect()
        return None

    # Iteratively filter beads with fewer than n_neighbors
    def filter_toofewneighbors(self):
        print(f"\nFiltering beads with <{self.n_neighbors} neighbors...")
        while any(np.diff(self.matrix.indptr) < self.n_neighbors):
            m = np.diff(self.matrix.indptr) >= self.n_neighbors
            self.apply_mask(m, "Too few neighbors")

    # Filter beads disconnected from the main graph
    def filter_disconnected(self):
        print(f"\nRemoving disconnected components...")
        n_components, labels = sp.csgraph.connected_components(self.matrix_k(), directed=False, return_labels=True)
        m = labels == np.bincount(labels).argmax()
        if np.sum(~m) > 0:
            self.apply_mask(m, "Disconnected")
    
    # Convert knn_matrix -> (knn_indices, knn_dists)
    def inddist(self, K=None):
        self.check()
        if K is None:
            K = self.n_neighbors
        
        knn_indices = np.zeros((self.matrix.shape[0], K), dtype=int)
        knn_dists = np.zeros((self.matrix.shape[0], K), dtype=float)
        
        for row_id in range(self.matrix.shape[0]):
            row_data = self.matrix[row_id].data
            row_indices = self.matrix[row_id].indices
            assert len(row_data) >= K
            row_nn_data_indices = np.argsort(row_data)[: K]
            knn_indices[row_id] = row_indices[row_nn_data_indices]
            knn_dists[row_id] = row_data[row_nn_data_indices]

        # check output?
        return knn_indices, knn_dists

    # Return the matrix subsetted to n_neighbors    
    def matrix_k(self):
        self.check()
        rows, cols, data = [], [], []
        for i in range(self.matrix.shape[0]):
            row_start, row_end = self.matrix.indptr[i], self.matrix.indptr[i+1]
            row_data = self.matrix.data[row_start:row_end]
            row_indices = self.matrix.indices[row_start:row_end]
    
            if len(row_data) > self.n_neighbors:
                top_k_indices = np.argpartition(row_data, self.n_neighbors)[:self.n_neighbors]
                row_data = row_data[top_k_indices]
                row_indices = row_indices[top_k_indices]
    
            rows.extend([i] * len(row_data))
            cols.extend(row_indices)
            data.extend(row_data)

        # check output? sort output?
        return sp.csr_matrix((data, (rows, cols)), shape=self.matrix.shape, dtype=np.float32)

    # Return an undirected igraph, weighted by cosine similarity
    def graph(self):
        import igraph as ig
        matrix_k = self.matrix_k()
        G = ig.Graph(n = matrix_k.shape[0],
                     edges = zip(*matrix_k.nonzero()),
                     edge_attrs = {'weight': 1-matrix_k.data}, # cosine similarity weights
                     directed = False)
        G.simplify(multiple=True, loops=True, combine_edges="max")
        return G

### LEIDEN METHODS #############################################################

def leiden_init(G, K=15, cores=1):
    print(f"    Leiden cluster...")
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
    import leidenalg as la
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
    print(f"    Leiden graph...")
    G.contract_vertices(membership, combine_attrs=None)
    G.simplify(multiple=True, loops=True, combine_edges=sum)
    
    # Convert Leiden graph to symmetric CSR
    i, j = zip(*G.get_edgelist())
    ic_edges = sp.csr_matrix((G.es["weight"], (i,j)), shape=(G.vcount(), G.vcount()))
    ic_edges = ic_edges.maximum(ic_edges.T)
    assert ic_edges.has_canonical_format
    assert np.all(ic_edges.diagonal() == 0)
    ic_edges.eliminate_zeros()

    # Create the Leiden KNN
    print(f"    Leiden KNN...")
    def my_knn(mat):
        from sklearn.preprocessing import normalize
        from sparse_dot_topn import sp_matmul_topn
        mat_norm = normalize(mat, norm='l2', axis=1, copy=True)
        C = sp_matmul_topn(mat_norm, mat_norm.T, top_n=K+1, n_threads=cores)
        C.data = 1 - C.data # cosine similarity -> distance
        C.data = np.maximum(C.data, 0)
        C.data = np.minimum(C.data, 1)
        C.setdiag(0)
        C.eliminate_zeros()
        # assert isinstance(C, sp.csr_matrix) and np.all(np.diff(C.indptr) == 14)
        return(C)
    
    leiden_knn_matrix = my_knn(ic_edges)

    # Create knn_indices, knn_dists
    lil = leiden_knn_matrix.tolil(copy=True)
    knn_indices = np.full((lil.shape[0], K), -1, dtype=np.int32)
    knn_dists = np.full((lil.shape[0], K), np.inf, dtype=np.float32)
    for i in range(lil.shape[0]):
        inds = np.array(lil.rows[i], dtype=np.int32)
        vals = np.array(lil.data[i], dtype=np.float32)
        sorted_indices = np.argsort(vals)
        inds = inds[sorted_indices][:K]
        vals = vals[sorted_indices][:K]
        knn_indices[i, :len(inds)] = inds
        knn_dists[i, :len(vals)] = vals
    
    print(f"    Leiden UMAP...")
    from umap import UMAP
    reducer = UMAP(n_components = 2,
                   metric = "cosine",
                   precomputed_knn = (knn_indices, knn_dists),
                   random_state = None,
                   verbose = False,
                   init = "spectral",
                   n_jobs = cores,
                   n_epochs = 5000,
                   n_neighbors=K,
                   min_dist=0.1
                  )
    embedding = reducer.fit_transform(ic_edges)
    assert np.all(np.isfinite(embedding))
    
    embedding -= np.mean(embedding, axis=0)
    return membership, ic_edges, embedding


def save_leiden_plots(embedding, knn, ic_edges, out_dir):
    print("\nLeiden plot 1...")
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.scatter(x=embedding[:, 0], y=embedding[:, 1], s=1)
    ax.set_title('Leiden Initialization')
    ax.axis('equal')
    fig.savefig(os.path.join(out_dir, "init.pdf"), dpi=200)

    print("\nLeiden plot 2...")
    fig, axs = plt.subplots(3, 2, figsize=(8, 8))
    def my_plot(embedding, vec, lab, axs):
        axs[0].hist(vec, bins=30, color='skyblue', edgecolor='black')
        axs[0].set_title(lab)
    
        axs[1].scatter(x=embedding[:,0], y=embedding[:,1], c=vec, cmap='viridis', s=1)
        axs[1].set(xticks=[], yticks=[], xticklabels=[], yticklabels=[])
        [spine.set_visible(False) for spine in axs[1].spines.values()]
        axs[1].set_aspect('equal')
        axs[1].set_title(lab)
    
    vec_beads = np.bincount(knn.membership)
    vec_strength = np.sum(ic_edges, axis=0).A1
    vec_logumi = np.log10(np.bincount(knn.membership, weights=knn.umi))
    
    my_plot(embedding, vec_beads, 'Beads', axs[0])
    my_plot(embedding, vec_strength, 'Strength', axs[1])
    my_plot(embedding, vec_logumi, 'log10(UMI)', axs[2])
    
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "init2.pdf"), dpi=200)


def save_umap_neighbor_plots(embedding, knn, out_dir):
    knn_indices, knn_dists = knn.inddist()
    
    # 1) Plot UMAP embedding neighbor distances
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    def my_hist(ax, dists, nn):
        ax.hist(np.log10(dists[:,nn][dists[:,nn] > 0]), bins=100)
        ax.set(xlabel='Distance (log10)', ylabel='Count', title=f'Distance to neighbor {nn+1}')
        meanval = np.mean(dists[:,nn])
        ax.axvline(np.log10(meanval), color='red', linestyle='dashed')
        ax.text(np.log10(meanval)+0.1, ax.get_ylim()[1] * 0.95, f'Mean: {meanval:.2f}', color='black', ha='left')

    dists = np.linalg.norm(embedding[knn_indices] - embedding[:,None,:], axis=2)
    my_hist(axes[0,0], dists, 0)
    my_hist(axes[0,1], dists, dists.shape[1]-1)
    my_hist(axes[1,0], dists, dists.shape[1]//2)
    
    axes[1,1].hexbin(x=np.log10(dists[dists > 0]), y=knn_dists[dists > 0],
                     gridsize=100, bins='log', cmap='plasma')
    axes[1,1].set_xlabel('UMAP Embedding distance (log10)')
    axes[1,1].set_ylabel('Cosine Distance')
    axes[1,1].set_title(f'Cosine vs. UMAP Distance ({knn.n_neighbors})')
    
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "umap_dists.pdf"), dpi=200)

    # 2) Plot UMAP embedding neighborhoods
    fig, axs = plt.subplots(4, 4, figsize=(8, 8))
    axs = axs.flatten() if type(axs) == np.ndarray else np.array([axs])
    for i,ind in enumerate(np.random.choice(embedding.shape[0], size=16, replace=False)):
        neighbors = embedding[knn_indices[ind]]
        allbeads = embedding[(embedding[:,0] >= np.min(neighbors[:,0])) & 
                             (embedding[:,0] <= np.max(neighbors[:,0])) &
                             (embedding[:,1] >= np.min(neighbors[:,1])) &
                             (embedding[:,1] <= np.max(neighbors[:,1]))]
        if allbeads.shape[0] > 100_000/16:
            # Don't plot neighborhoods that are excessively large (too many background points)
            axs[i].text(x=0.5, y=0.5, s='Neighborhood\ntoo large',
                horizontalalignment='center', verticalalignment='center',
                transform=axs[i].transAxes, fontsize=10, color='red')
            axs[i].set_axis_off()
        else:
            # Plot the neighborhood of a random bead (red) colored by cosine similarity
            axs[i].scatter(allbeads[:,0], allbeads[:,1], color='grey', s=15, alpha=0.5)
            axs[i].scatter(neighbors[:,0], neighbors[:,1], c=1-knn_dists[ind,:], s=15, cmap='viridis', vmin=0, vmax=1)
            axs[i].scatter(embedding[ind,0], embedding[ind,1], color='red', s=15)
            axs[i].set_aspect('equal', adjustable='box')
            axs[i].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
            axs[i].set_title(knn.sb[ind], fontsize=8)
    
    fig.savefig(os.path.join(out_dir, "umap_neighbors.pdf"), dpi=200)

