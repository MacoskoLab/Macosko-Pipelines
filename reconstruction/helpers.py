import os
import math
import heapq
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
from functools import reduce
from collections import Counter

# (x, y)
def hexmap(embedding, title="", fontsize=12):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    x, y = embedding[:, 0], embedding[:, 1]
    hb = ax.hexbin(x, y, cmap='viridis', linewidths=0.5)
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
    n = math.ceil(len(embeddings)**0.5)
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

# (sb1, sb2, umi)
def uvc(df):
    assert df.shape[1] == 3
    df.columns = ['sb1_index', 'sb2_index', 'umi']
    
    df1 = df.groupby('sb1_index').agg(umi=('umi', 'sum'), size=('sb1_index', 'size')).reset_index()
    df1[['umi', 'size']] = np.log10(df1[['umi', 'size']])
    
    df2 = df.groupby('sb2_index').agg(umi=('umi', 'sum'), size=('sb2_index', 'size')).reset_index()
    df2[['umi', 'size']] = np.log10(df2[['umi', 'size']])

    def plot(ax, df, title):
        hb = ax.hexbin(df['umi'], df['size'], gridsize=100, bins='log', cmap='plasma')
        ax.set_xlabel(f'umi (mean: {np.mean(df["umi"]):.2f}, median: {np.median(df["umi"]):.2f})')
        ax.set_ylabel(f'connections (mean: {np.mean(df["size"]):.2f}, median: {np.median(df["size"]):.2f})')
        ax.set_title(title)

        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        pts = np.linspace(max(x1,y1), min(x2,y2), 100)
        ax.plot(pts, pts, color='black', linewidth=0.5)
        ax.set_xlim((x1, x2))
        ax.set_ylim((y1, y2))
        ax.axis('equal')

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    plot(axes[0], df1, "sb1")
    plot(axes[1], df2, "sb2")
    plt.tight_layout()
    
    return fig, axes

# (x, y, color)
def beadplot(puck, cmap='viridis'):
    if isinstance(puck, np.ndarray):
        if puck.shape[1] == 2:
            puck = np.hstack((puck, np.zeros((puck.shape[0], 1))))
        puck = puck[puck[:, 2].argsort()]
        x = puck[:,0]
        y = puck[:,1]
        umi = puck[:,2]
    elif isinstance(puck, pd.DataFrame):
        if puck.shape[1] == 2:
            puck.loc[:,"color"] = 0
        puck = puck.sort_values(by=puck.columns[2])
        x = puck.iloc[:,0]
        y = puck.iloc[:,1]
        umi = puck.iloc[:,2]
    else:
        raise TypeError("Input must be a NumPy ndarray or a pandas DataFrame")

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, c=umi, cmap=cmap, s=0.1)
    plt.colorbar()
    plt.xlabel('xcoord')
    plt.ylabel('ycoord')
    plt.axis('square')
    return plt

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
    
### BEAD FILTERING METHODS #####################################################

def connection_filter(df):
    assert all(df.columns == ["sb1_index", "sb2_index", "umi"])
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    meta = {"umi_init": sum(df["umi"])}

    def hist_z(ax, data, z_high=None, z_low=None):
        ax.hist(data, bins=100)
        ax.set_ylabel('Count')
    
        meanval = np.mean(data)
        ax.axvline(meanval, color='black', linestyle='dashed')
        ax.text(meanval, ax.get_ylim()[1] * 0.95, f'Mean: {10**meanval:.2f}', color='black', ha='center')
    
        if z_high:
            assert z_high > 0
            zval = meanval + np.std(data) * z_high
            ax.axvline(zval, color='red', linestyle='dashed')
            ax.text(zval, ax.get_ylim()[1]*0.9, f'{z_high}z: {10**zval:.2f}', color='red', ha='left')
    
        if z_low:
            assert z_low < 0
            zval = meanval + np.std(data) * z_low
            ax.axvline(zval, color='red', linestyle='dashed')
            ax.text(zval, ax.get_ylim()[1]*0.9, f'{z_low}z: {10**zval:.2f}', color='red', ha='right')
    
    # Remove sb1 beads
    z_high = 3
    sb1 = df.groupby(f'sb1_index').agg(umi=('umi', 'sum'), connections=('umi', 'size'), max=('umi', 'max')).reset_index()
    hist_z(axes[0,0], np.log10(sb1['connections']), z_high)
    axes[0,0].set_xlabel('Connections')
    axes[0,0].set_title('sb1 connections')
    hist_z(axes[0,1], np.log10(sb1['max']))
    axes[0,1].set_xlabel('Max')
    axes[0,1].set_title('sb1 max')

    logcon = np.log10(sb1['connections'])
    high = np.where(logcon >= np.mean(logcon) + np.std(logcon) * z_high)[0]
    sb1_remove = reduce(np.union1d, [high])
    df = df[~df['sb1_index'].isin(sb1_remove)]
    
    meta["sb1_high"] = len(high)
    meta["sb1_removed"] = len(sb1_remove)
    meta["umi_half"] = sum(df["umi"])
    print(f"{len(high)} high R1 beads ({len(high)/len(sb1)*100:.2f}%)")
    diff = meta['umi_init']-meta['umi_half'] ; print(f"{diff} R1 UMIs filtered ({diff/meta['umi_init']*100:.2f}%)")

    # Remove sb2 beads
    z_high = 3 #; z_low = -3
    sb2 = df.groupby(f'sb2_index').agg(umi=('umi', 'sum'), connections=('umi', 'size'), max=('umi', 'max')).reset_index()
    hist_z(axes[1,0], np.log10(sb2['connections']), z_high) #, z_low)
    axes[1,0].set_xlabel('Connections')
    axes[1,0].set_title('sb2 connections')
    hist_z(axes[1,1], np.log10(sb2['max']))
    axes[1,1].set_xlabel('Max')
    axes[1,1].set_title('sb2 max')
    
    logcon = np.log10(sb2['connections'])
    high = np.where(logcon >= np.mean(logcon) + np.std(logcon) * z_high)[0]
    # low = np.where(logcon <= np.mean(logcon) + np.std(logcon) * z_low)[0]
    # noise = np.where(sb2['max'] <= 1)[0]
    sb2_remove = reduce(np.union1d, [high]) # [high, low, noise]
    df = df[~df['sb2_index'].isin(sb2_remove)]
    
    meta["sb2_high"] = len(high)
    # meta["sb2_low"] = len(low)
    # meta["sb2_noise"] = len(noise)
    meta["sb2_removed"] = len(sb2_remove)
    meta["umi_final"] = sum(df["umi"])
    print(f"{len(high)} high R2 beads ({len(high)/len(sb2)*100:.2f}%)")
    # print(f"{len(low)} low R2 beads ({len(low)/len(sb2)*100:.2f}%)")
    # print(f"{len(noise)} noise R2 beads ({len(noise)/len(sb2)*100:.2f}%)")
    diff = meta['umi_half']-meta['umi_final'] ; print(f"{diff} R2 UMIs filtered ({diff/meta['umi_half']*100:.2f}%)")
    
    # Factorize the new dataframe
    codes1, uniques1 = pd.factorize(df['sb1_index'], sort=True)
    df.loc[:, 'sb1_index'] = codes1
    codes2, uniques2 = pd.factorize(df['sb2_index'], sort=True)
    df.loc[:, 'sb2_index'] = codes2

    # assert set(df.sb1_index) == set(range(max(df.sb1_index)+1))
    # assert set(df.sb2_index) == set(range(max(df.sb2_index)+1))
    
    fig.tight_layout()
    return df, uniques1, uniques2, fig, meta


def knn_filter(knn_indices, knn_dists):
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    filter_indexes = set()
    meta = dict()

    def hist_z(ax, data, z_high=None, z_low=None, bins=100):
        ax.hist(data, bins=bins)
        ax.set_ylabel('Count')
    
        meanval = np.mean(data)
        ax.axvline(meanval, color='black', linestyle='dashed')
        ax.text(meanval, ax.get_ylim()[1] * 0.95, f'Mean: {meanval:.2f}', color='black', ha='center')
    
        if z_high:
            assert z_high > 0
            zval = meanval + np.std(data) * z_high
            ax.axvline(zval, color='red', linestyle='dashed')
            ax.text(zval, ax.get_ylim()[1]*0.9, f'{z_high}z: {zval:.2f}', color='red', ha='left')
    
        if z_low:
            assert z_low < 0
            zval = meanval + np.std(data) * z_low
            ax.axvline(zval, color='red', linestyle='dashed')
            ax.text(zval, ax.get_ylim()[1]*0.9, f'{z_low}z: {zval:.2f}', color='red', ha='right')
    
    # Filter beads with far nearest neighbors
    z_high = 3
    data = knn_dists[:,1]
    hist_z(axes[0,0], data, z_high)
    axes[0,0].set_xlabel(f'Distance')
    axes[0,0].set_title(f'Nearest neighbor distance')

    high = knn_indices[data >= np.mean(data) + np.std(data) * z_high, 0]
    filter_indexes.update(high)
    print(f"{len(high)} far-NN beads removed")
    meta["far-NN"] = len(high)

    # Plot furthest neighbor distance
    hist_z(axes[0,1], knn_dists[:,-1])
    axes[0,1].set_xlabel(f'Distance')
    axes[0,1].set_title(f'Furthest neighbor distance ({knn_dists.shape[1]})')
    
    # Filter too-high or too-low in-edges
    z_high = 3 ; z_low = -3
    indexes, data = np.unique(knn_indices, return_counts=True)
    hist_z(axes[1,0], data, z_high, z_low, bins=np.arange(0, data.max()+1))
    axes[1,0].set_xlabel('In-edges')
    axes[1,0].set_title('Number of in-edges')
    axes[1,0].set_xlim(0, 4*knn_indices.shape[1])
    
    high = indexes[data >= np.mean(data) + np.std(data) * z_high]
    filter_indexes.update(high)
    print(f"{len(high)} in-high beads removed")
    meta["in-high"] = len(high)
    
    low = indexes[data <= np.mean(data) + np.std(data) * z_low]
    filter_indexes.update(low)
    print(f"{len(low)} in-low beads removed")
    meta["in-low"] = len(low)

    # Filter low clustering coefficients
    z_low = -3
    knn_matrix = create_knn_matrix(knn_indices, knn_dists)
    G = nx.from_scipy_sparse_array(knn_matrix, create_using=nx.Graph, edge_attribute=None) # undirected, unweighted
    clustering = nx.clustering(G, nodes=None, weight=None)
    data = [clustering[key] for key in sorted(clustering.keys())]
    hist_z(axes[1,1], data, z_low=z_low)
    axes[1,1].set_xlabel('Clustering coefficient')
    axes[1,1].set_title('Local clustering coefficient')
    
    low = knn_indices[data <= np.mean(data) + np.std(data) * z_low, 0]
    filter_indexes.update(low)
    print(f"{len(low)} cluster-low beads removed")
    meta["cluster-low"] = len(low)

    # Filter weakly-connected components
    n_components, labels = sp.csgraph.connected_components(csgraph=knn_matrix, directed=True, connection='strong')
    wcc = np.where(labels != np.bincount(labels).argmax())[0]
    filter_indexes.update(wcc)
    print(f"{len(wcc)} WCC beads removed")
    meta["wcc"] = len(wcc)
     
    fig.tight_layout()
    return filter_indexes, fig, meta

### KNN METHODS ################################################################

# Compute the KNN using NNDescent
def knn_descent(mat, n_neighbors, metric="cosine", n_jobs=-1):    
    from pynndescent import NNDescent
    knn_search_index = NNDescent(
        data=mat,
        metric=metric,
        metric_kwds={},
        n_neighbors=n_neighbors,
        n_trees=64, # originally None
        # leaf_size=None,
        pruning_degree_multiplier=3.0, # originally 1.5
        diversify_prob=0.0, # originally 1.0
        # tree_init=True,
        # init_graph=None,
        # init_dist=None,
        random_state=None,
        low_memory=True, # originally False
        max_candidates=60, # originally None
        max_rptree_depth=999999, # originally 200
        n_iters=512, # originally None
        delta=0.0001, # originally 0.001
        n_jobs=n_jobs,
        # compressed=False,
        # parallel_batch_queries=False,
        verbose=True # originally False
    )
    knn_indices, knn_dists = knn_search_index.neighbor_graph
    return knn_indices, knn_dists

# UMAP NNDescent wrapper
def nearest_neighbors(mat, n_neighbors, metric="cosine", n_jobs=-1):
    from umap.umap_ import nearest_neighbors
    knn_indices, knn_dists, _ = nearest_neighbors(mat,
                                    n_neighbors = n_neighbors,
                                    metric = metric,
                                    metric_kwds = {},
                                    angular = False, # Does nothing?
                                    random_state = None,
                                    low_memory = True, # False?
                                    use_pynndescent = True, # Does nothing?
                                    n_jobs = n_jobs,
                                    verbose = True
                                )
    return knn_indices, knn_dists

# knn_indices1 can be row-sliced, knn_indices2 is original
def knn_merge(knn_indices1, knn_dists1, knn_indices2, knn_dists2):
    assert knn_indices1.shape == knn_dists1.shape
    assert knn_indices2.shape == knn_dists2.shape
    assert knn_indices1.shape[0] == knn_dists1.shape[0] == knn_indices2.shape[0] == knn_dists2.shape[0]

    assert all(knn_indices2[:,0] == np.arange(knn_indices2.shape[0]))
    assert np.all(knn_indices2 >= 0) and np.all(knn_dists2 >= 0)
    index_map = dict(zip(knn_indices1[:,0], knn_indices2[:,0]))
    knn_indices1 = np.vectorize(lambda i: index_map.get(i,-1))(knn_indices1)
    assert all(knn_indices1[:,0] == knn_indices2[:,0])
    # assert np.all(knn_dists1[np.where(knn_indices1==knn_indices2)] == knn_dists2[np.where(knn_indices1==knn_indices2)])
    
    k = max(knn_indices1.shape[1], knn_indices2.shape[1])
    knn_indices = np.zeros((knn_indices1.shape[0], k), dtype=np.int32)
    knn_dists = np.zeros((knn_indices1.shape[0], k), dtype=np.float64)
    
    for i in range(knn_indices1.shape[0]):
        d1 = {i:d for i,d in zip(knn_indices1[i], knn_dists1[i]) if i >= 0}
        d2 = {i:d for i,d in zip(knn_indices2[i], knn_dists2[i]) if i >= 0}
        d = d1 | d2

        row = sorted(d.items(), key=lambda item: item[1])[:k]
        inds, dists = zip(*row)
        
        knn_indices[i,:] = inds
        knn_dists[i,:] = dists

    return knn_indices, knn_dists

### MNN METHODS ################################################################
### source: https://umap-learn.readthedocs.io/en/latest/mutual_nn_umap.html ####

def create_knn_matrix(knn_indices, knn_dists):
    assert knn_indices.shape == knn_dists.shape
    assert knn_indices.dtype == np.int32 and knn_dists.dtype == np.float64
    assert np.max(knn_indices) < len(knn_indices)
    assert np.array_equal(knn_indices[:,0], np.arange(len(knn_indices)))
    assert np.all(knn_dists[:,1:] > 0) and not np.any(np.isnan(knn_dists))
    
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

def min_spanning_tree(knn_matrix):
    Tcsr = sp.csgraph.minimum_spanning_tree(knn_matrix)
    Tcsr = sp.coo_matrix(Tcsr)
    weights_tuples = zip(Tcsr.row, Tcsr.col, Tcsr.data)
    sorted_weights_tuples = sorted(weights_tuples, key=lambda tup: tup[2])
    return sorted_weights_tuples 

# Search to find path neighbors (todo: stop at dist=0)
def find_new_nn(indices, dists, out_neighbors, i_range):
    mnn_dists = [] 
    mnn_indices = []
    
    for i in i_range:
        min_indices = []
        min_distances = []
        
        heap = [(0,i)] ; heapq.heapify(heap) 
        mapping = {}
        seen = set()
        
        while len(min_distances) < out_neighbors and len(heap) > 0:
            dist, nn = heapq.heappop(heap)
            if nn in seen or nn < 0:
                continue
        
            min_indices.append(nn)
            min_distances.append(dist)
            seen.add(nn)
            
            for nn_i, nn_d in zip(indices[nn], dists[nn]):
                if nn_i in seen or nn_d <= 0:
                    continue
                distance = dist + nn_d
                if nn_i not in mapping or distance < mapping[nn_i]:
                    mapping[nn_i] = distance
                    heapq.heappush(heap, (distance, nn_i))
            
        if len(min_distances) < out_neighbors:
            for i in range(out_neighbors-len(min_distances)):
                min_indices.append(-1)
                min_distances.append(np.inf)
        
        mnn_indices.append(min_indices)
        mnn_dists.append(min_distances)
        
    return np.array(mnn_indices, dtype=np.int32), np.array(mnn_dists)

# Prune non-reciprocated edges and add MST edges
def create_mnn(knn_indices, knn_dists, in_neighbors):
    # Create mutual graph
    print("Creating mutual graph...")
    knn_matrix = create_knn_matrix(knn_indices, knn_dists).tocsr()
    lil = knn_matrix.multiply(knn_matrix.transpose()).sqrt().tolil()
    lil.setdiag(0)
    
    # Add MST edges
    print("Adding MST edges...")
    sorted_weights_tuples = min_spanning_tree(knn_matrix)
    for i,j,v in sorted_weights_tuples:
        lil[i,j] = v
        lil[j,i] = v
    del sorted_weights_tuples, knn_matrix
    
    # Memory-map the data
    print("Memory-mapping the data...")
    indices = np.memmap('/dev/shm/indices.bin', dtype='int32', mode='w+', shape=(knn_indices.shape[0],in_neighbors))
    dists = np.memmap('/dev/shm/dists.bin', dtype='float64', mode='w+', shape=(knn_dists.shape[0],in_neighbors))
    indices[:] = -1 ; dists[:] = -1
    for i in range(lil.shape[0]):
        cols = np.array(lil.rows[i], dtype=np.int32)
        vals = np.array(lil.data[i], dtype=np.float64)
        if len(cols) > in_neighbors:
            inds = np.argpartition(vals, in_neighbors)[:in_neighbors]
            cols = cols[inds]
            vals = vals[inds]

        indices[i, :len(cols)] = cols
        dists[i, :len(vals)] = vals
    
    indices.flush() ; dists.flush()
    import gc ; gc.collect()
    del indices, dists
    import gc ; gc.collect()
    shape = (knn_indices.shape[0], in_neighbors)
    return shape
    
def find_path_neighbors(shape, out_neighbors, n_jobs=-1):
    print("Finding new path neighbors...")
    from multiprocessing import Pool
    if n_jobs < 1:
        n_jobs = len(os.sched_getaffinity(0))
        
    indices = np.memmap('/dev/shm/indices.bin', dtype='int32', mode='r', shape=shape)
    dists = np.memmap('/dev/shm/dists.bin', dtype='float64', mode='r', shape=shape)
    ranges = np.array_split(range(len(indices)), n_jobs)
    with Pool(processes=n_jobs) as pool:
        results = pool.starmap(find_new_nn, [(indices, dists, out_neighbors, i_range) for i_range in ranges])

    mnn_indices, mnn_dists = zip(*results)
    mnn_indices = np.vstack(mnn_indices)
    mnn_dists = np.vstack(mnn_dists)
    
    assert mnn_indices.shape == mnn_dists.shape
    assert mnn_indices.shape[0] == shape[0]
    assert np.max(mnn_indices) < len(mnn_indices)
    assert np.all(mnn_indices >= 0) and np.all(mnn_dists >= 0)
    assert np.array_equal(mnn_indices[:,0], np.arange(len(mnn_indices)))
    assert mnn_indices.dtype == np.int32 and mnn_dists.dtype == np.float64
    return mnn_indices, mnn_dists
