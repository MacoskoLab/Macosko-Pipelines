import math
import numpy as np
import pandas as pd
import networkx as nx
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
    n = math.ceil(len(embeddings)**0.5)
    fig, axes = plt.subplots(n, n, figsize=(8, 8))
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

def procrustes_distance(p1, p2):
    from scipy.spatial import procrustes
    assert p1.shape == p2.shape
    _, _, disparity = procrustes(p1, p2)
    return disparity

def convergence_plot(embeddings):
    assert len(embeddings) > 1
    x = [(i+2)*1000 for i in range(len(embeddings)-1)]
    y1 = [L2_distance(embeddings[i], embeddings[i+1]) for i in range(len(embeddings)-1)]
    y2 = [procrustes_distance(embeddings[i], embeddings[i+1]) for i in range(len(embeddings)-1)]

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    axes[0].scatter(x, y1, color='blue')
    axes[0].plot(x, y1, color='red')
    axes[0].set_xlabel('Epochs')
    axes[0].set_ylabel('Mean UMAP Distance')
    axes[0].set_title('Average L2 Displacement')
    
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
    z_high = 3 ; z_low = -3
    sb2 = df.groupby(f'sb2_index').agg(umi=('umi', 'sum'), connections=('umi', 'size'), max=('umi', 'max')).reset_index()
    hist_z(axes[1,0], np.log10(sb2['connections']), z_high, z_low)
    axes[1,0].set_xlabel('Connections')
    axes[1,0].set_title('sb2 connections')
    hist_z(axes[1,1], np.log10(sb2['max']))
    axes[1,1].set_xlabel('Max')
    axes[1,1].set_title('sb2 max')
    
    logcon = np.log10(sb2['connections'])
    high = np.where(logcon >= np.mean(logcon) + np.std(logcon) * z_high)[0]
    low = np.where(logcon <= np.mean(logcon) + np.std(logcon) * z_low)[0]
    noise = np.where(sb2['max'] <= 1)[0]
    sb2_remove = reduce(np.union1d, [high, low, noise])
    df = df[~df['sb2_index'].isin(sb2_remove)]
    
    meta["sb2_high"] = len(high)
    meta["sb2_low"] = len(low)
    meta["sb2_noise"] = len(noise)
    meta["sb2_removed"] = len(sb2_remove)
    meta["umi_final"] = sum(df["umi"])
    print(f"{len(high)} high R2 beads ({len(high)/len(sb2)*100:.2f}%)")
    print(f"{len(low)} low R2 beads ({len(low)/len(sb2)*100:.2f}%)")
    print(f"{len(noise)} noise R2 beads ({len(noise)/len(sb2)*100:.2f}%)")
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
    knn_matrix = create_knn_matrix(knn_indices, knn_dists, knn_indices.shape[1])
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
    n_components, labels = scipy.sparse.csgraph.connected_components(csgraph=knn_matrix, directed=True, connection='strong')
    wcc = np.where(labels != np.bincount(labels).argmax())[0]
    filter_indexes.update(wcc)
    print(f"{len(wcc)} WCC beads removed")
    meta["wcc"] = len(wcc)
     
    fig.tight_layout()
    return filter_indexes, fig, meta

### KNN METHODS ################################################################

def knn_descent(mat, n_neighbors, metric="cosine", n_cores=-1):
    # from umap.umap_ import nearest_neighbors
    # knn_indices, knn_dists, _ = nearest_neighbors(mat,
    #                                 n_neighbors = n_neighbors,
    #                                 metric = metric,
    #                                 metric_kwds = {},
    #                                 angular = False, # Does nothing?
    #                                 random_state = None,
    #                                 low_memory = True, # False?
    #                                 use_pynndescent = True, # Does nothing?
    #                                 n_jobs = n_cores,
    #                                 verbose = True
    #                             )
    
    from pynndescent import NNDescent
    knn_search_index = NNDescent(
        data=mat,
        metric=metric,
        metric_kwds={},
        n_neighbors=n_neighbors,
        n_trees=64, # use theirs? nah bigger?
        leaf_size=None, # bigger?
        pruning_degree_multiplier=1.5, # 3.0
        diversify_prob=1.0, # 0.0
        # tree_init=True,
        # init_graph=None,
        # init_dist=None,
        random_state=None,
        low_memory=True, # False?
        max_candidates=None, # no touchy
        max_rptree_depth=200, # 100?
        n_iters=100, # make bigger? use theirs?
        delta=0.001, # don't touch? or make smaller?
        n_jobs=-1,
        compressed=True, # True?
        parallel_batch_queries=False, # True?
        verbose=True
    )
    knn_indices, knn_dists = knn_search_index.neighbor_graph

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
import scipy

def create_knn_matrix(knn_indices, knn_dists, n_neighbors):
    assert knn_indices.shape == knn_dists.shape
    assert n_neighbors <= knn_indices.shape[1]
    rows = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
    cols = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
    vals = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.float32)
    
    pos = 0
    for i, indices in enumerate(knn_indices):
        for j, index in enumerate(indices[:n_neighbors]):
            if index == -1:
                continue
            rows[pos] = i 
            cols[pos] = index
            vals[pos] = knn_dists[i][j]
            pos += 1
    
    knn_matrix = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
    return knn_matrix
    
def min_spanning_tree(knn_matrix):
    Tcsr = scipy.sparse.csgraph.minimum_spanning_tree(knn_matrix)
    Tcsr = scipy.sparse.coo_matrix(Tcsr)
    weights_tuples = zip(Tcsr.row, Tcsr.col, Tcsr.data)
    sorted_weights_tuples = sorted(weights_tuples, key=lambda tup: tup[2])
    return sorted_weights_tuples 

def create_connected_graph(mutual_nn, total_mutual_nn, knn_indices, knn_dists, n_neighbors, connectivity):
    import copy
    connected_mnn = copy.deepcopy(mutual_nn)
    assert connectivity in ["min_tree", "full_tree"] # "nearest" produced np.inf
    
    if connectivity == "nearest":
        for i in range(len(knn_indices)): 
            if len(mutual_nn[i]) == 0:
                first_nn = knn_indices[i][1]
                if first_nn != -1:
                    connected_mnn[i].add(first_nn) 
                    connected_mnn[first_nn].add(i) 
                    total_mutual_nn += 1
        return connected_mnn
    
    # Create graph for mutual NN
    rows = np.zeros(total_mutual_nn, dtype=np.int32)
    cols = np.zeros(total_mutual_nn, dtype=np.int32)
    vals = np.zeros(total_mutual_nn, dtype=np.float32)
    pos = 0
    for i in connected_mnn:
        for j in connected_mnn[i]:
            rows[pos] = i 
            cols[pos] = j
            vals[pos] = 1
            pos += 1
    graph = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
    
    # Find number of connected components
    n_components, labels = scipy.sparse.csgraph.connected_components(csgraph=graph, directed=True, return_labels=True, connection='strong')
    print(f"connected_components: {n_components}")
    label_mapping = {i:[] for i in range(n_components)}
    
    for index, component in enumerate(labels):
        label_mapping[component].append(index)
    
    # Find the min spanning tree with KNN
    sorted_weights_tuples = min_spanning_tree(create_knn_matrix(knn_indices, knn_dists, n_neighbors))
    
    # Add edges until graph is connected
    for pos,(i,j,v) in enumerate(sorted_weights_tuples):
        
        if connectivity == "full_tree":
            connected_mnn[i].add(j)
            connected_mnn[j].add(i) 
          
        elif connectivity == "min_tree" and labels[i] != labels[j]:
            if len(label_mapping[labels[i]]) < len(label_mapping[labels[j]]):
                i, j = j, i
              
            connected_mnn[i].add(j)
            connected_mnn[j].add(i)
            j_pos = label_mapping[labels[j]]
            labels[j_pos] = labels[i]
            label_mapping[labels[i]].extend(j_pos)
    
    return connected_mnn  

# Search to find path neighbors
def find_new_nn(knn_indices, knn_dists, knn_indices_pos, connected_mnn, n_neighbors_max):
    import heapq
    new_knn_dists = [] 
    new_knn_indices = []
    
    for i in range(len(knn_indices)): 
        min_distances = []
        min_indices = []
        
        heap = [(0,i)]
        mapping = {}
              
        seen = set()
        heapq.heapify(heap) 
        while(len(min_distances) < n_neighbors_max and len(heap) >0):
            dist, nn = heapq.heappop(heap)
            if nn == -1:
                continue
            
            if nn not in seen:
                min_distances.append(dist)
                min_indices.append(nn)
                seen.add(nn)
                neighbor = connected_mnn[nn]
                
                for nn_nn in neighbor:
                    if nn_nn not in seen:
                        distance = 0
                        if nn_nn in knn_indices_pos[nn]:
                            pos = knn_indices_pos[nn][nn_nn]
                            distance = knn_dists[nn][pos] 
                        else:
                            pos = knn_indices_pos[nn_nn][nn]
                            distance = knn_dists[nn_nn][pos] 
                        distance += dist
                        if nn_nn not in mapping:
                            mapping[nn_nn] = distance
                            heapq.heappush(heap, (distance, nn_nn))
                        elif mapping[nn_nn] > distance:
                            mapping[nn_nn] = distance
                            heapq.heappush(heap, (distance, nn_nn))
            
        if len(min_distances) < n_neighbors_max:
            for i in range(n_neighbors_max-len(min_distances)):
                min_indices.append(-1)
                min_distances.append(np.inf)
        
        new_knn_dists.append(min_distances)
        new_knn_indices.append(min_indices)
        
        if i % int(len(knn_dists) / 10) == 0:
            print("\tcompleted ", i, " / ", len(knn_dists), "epochs")
    return new_knn_dists, new_knn_indices

# Calculate the connected mutual nn graph
def mutual_nn_nearest(knn_indices, knn_dists, n_neighbors, n_neighbors_max, connectivity):
    mutual_nn = {}
    nearest_n = {}
    
    knn_indices_pos = [None] * len(knn_indices)
    for i, top_vals in enumerate(knn_indices):
        nearest_n[i] = set(top_vals)
        knn_indices_pos[i] = {}
        for pos, nn in enumerate(top_vals):
            knn_indices_pos[i][nn] = pos
    
    total_mutual_nn = 0
    for i, top_vals in enumerate(knn_indices):
        mutual_nn[i] = set()
        for ind, nn in enumerate(top_vals):
             if nn != -1 and (i in nearest_n[nn] and i != nn):
                 mutual_nn[i].add(nn)
                 total_mutual_nn += 1

    print("Creating connected graph...")
    connected_mnn = create_connected_graph(mutual_nn, total_mutual_nn, knn_indices, knn_dists, n_neighbors, connectivity)
    
    print("Finding new nearest neighbors...")
    new_knn_dists, new_knn_indices = find_new_nn(knn_indices, knn_dists, knn_indices_pos, connected_mnn, n_neighbors_max)
    
    return np.array(new_knn_indices, dtype=np.int32), np.array(new_knn_dists)
