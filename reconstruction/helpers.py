import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# (x, y)
def hexmap(embedding):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    x, y = embedding[:, 0], embedding[:, 1]
    hb = ax.hexbin(x, y, cmap='viridis')
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# [(x, y)]
def hexmaps(embeddings, titles = [], fontsize=10):
    n = math.ceil(len(embeddings)**0.5)
    fig, axes = plt.subplots(n, n, figsize=(2*n, 2*n))
    axes = axes.flatten()

    if len(titles) == 0:
        titles = ["" for _ in range(len(embeddings))]
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)
    assert len(titles) == len(embeddings)
    
    for ax, embedding, title in zip(axes, embeddings, titles):
        x, y = embedding[:, 0], embedding[:, 1]
        hb = ax.hexbin(x, y, cmap='viridis')
        ax.axis('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(f'{title}', fontsize=fontsize)
    for ax in axes[len(embeddings):]:
        ax.set_visible(False)

# (sb1, sb2, umi)
def uvc(df):
    assert df.shape[1] == 3
    df.columns = ['sb1', 'sb2', 'umi']
    
    df1 = df.groupby('sb1').agg(umi=('umi', 'sum'), size=('sb1', 'size')).reset_index()
    df1[['umi', 'size']] = np.log10(df1[['umi', 'size']])
    
    df2 = df.groupby('sb2').agg(umi=('umi', 'sum'), size=('sb2', 'size')).reset_index()
    df2[['umi', 'size']] = np.log10(df2[['umi', 'size']])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
    m = 1

    hb = ax1.hexbin(df1['umi'], df1['size'], gridsize=100, bins='log', cmap='plasma')
    ax1.set_xlabel(f'umi (mean: {np.mean(df1["umi"]):.2f}, median: {np.median(df1["umi"]):.2f})')
    ax1.set_ylabel(f'connections (mean: {np.mean(df1["size"]):.2f}, median: {np.median(df1["size"]):.2f})')
    ax1.set_title('sb1')

    xlim1 = ax1.get_xlim()
    ylim1 = ax1.get_ylim()
    x_vals1 = np.linspace(xlim1[0], xlim1[1], 100)
    ax1.plot(x_vals1, x_vals1 * m, color='black', linewidth=0.5, label=f'y = {m}x')
    ax1.set_xlim(xlim1)
    ax1.set_ylim(ylim1)
    
    hb = ax2.hexbin(df2['umi'], df2['size'], gridsize=100, bins='log', cmap='plasma')
    ax2.set_xlabel(f'umi (mean: {np.mean(df2["umi"]):.2f}, median: {np.median(df2["umi"]):.2f})')
    ax2.set_ylabel(f'connections (mean: {np.mean(df2["size"]):.2f}, median: {np.median(df2["size"]):.2f})')
    ax2.set_title('sb2')

    xlim2 = ax2.get_xlim()
    ylim2 = ax2.get_ylim()
    x_vals2 = np.linspace(xlim2[0], xlim2[1], 100)
    ax2.plot(x_vals2, x_vals2 * m, color='black', linewidth=0.5, label=f'y = {m}x')
    ax2.set_xlim(xlim2)
    ax2.set_ylim(ylim2)
    
    plt.tight_layout()
    plt.show()

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
        puck = puck.sort_values(by=puck.columns[2])
        x = puck.iloc[:,0]
        y = puck.iloc[:,1]
        umi = puck.iloc[:,2]
    else:
        raise TypeError("Input must be a NumPy ndarray or a pandas DataFrame")

    plt.figure(figsize=(12, 12))
    plt.scatter(x, y, c=umi, cmap=cmap, s=0.1)
    plt.colorbar()
    plt.xlabel('xcoord')
    plt.ylabel('ycoord')
    plt.axis('square')
    plt.show()

# embedding1, embedding2
def L2_distance(p1, p2):
    dists = np.linalg.norm(p1 - p2, axis=1)
    return np.sum(dists)

# embedding1, embedding2
def ICP_distance(points1, points2):
    import open3d as o3d
    points1_3d = np.hstack((points1, np.zeros((points1.shape[0], 1))))
    points2_3d = np.hstack((points2, np.zeros((points2.shape[0], 1))))
    
    pc1 = o3d.geometry.PointCloud()
    pc1.points = o3d.utility.Vector3dVector(points1_3d)
    pc2 = o3d.geometry.PointCloud()
    pc2.points = o3d.utility.Vector3dVector(points2_3d)
    
    threshold = 0.01
    trans_init = np.eye(4)
    transformation = o3d.pipelines.registration.registration_icp(
        pc1, pc2, threshold, trans_init,
        o3d.pipelines.registration.TransformationEstimationPointToPoint()
    ).transformation
    
    pc2.transform(transformation)
    distances = np.asarray(pc1.compute_point_cloud_distance(pc2))
    total_distance = np.sum(distances)
    return(total_distance)

### KNN METHODS ################################################################

def knn_descent(mat, n_neighbors, metric="cosine", n_cores=-1):
    from umap.umap_ import nearest_neighbors
    knn_indices, knn_dists, _ = nearest_neighbors(mat,
                                    n_neighbors = n_neighbors,
                                    metric = metric,
                                    metric_kwds = {},
                                    angular = False, # Does nothing?
                                    random_state = None, # sklearn.utils.check_random_state(0)
                                    low_memory = True, # False?
                                    use_pynndescent = True, # Does nothing?
                                    n_jobs = n_cores,
                                    verbose = True
                                )
    return knn_indices, knn_dists

### MNN METHODS ################################################################
### source: https://umap-learn.readthedocs.io/en/latest/mutual_nn_umap.html ####
import scipy

# Calculate min spanning tree
def min_spanning_tree(knn_indices, knn_dists, n_neighbors, threshold):
    rows = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
    cols = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
    vals = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.float32)
    
    pos = 0
    for i, indices in enumerate(knn_indices):
        for j, index in enumerate(indices[:threshold]):
            if index == -1:
                continue
            rows[pos] = i 
            cols[pos] = index
            vals[pos] = knn_dists[i][j]
            pos += 1
    
    matrix = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
    Tcsr = scipy.sparse.csgraph.minimum_spanning_tree(matrix)
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
    print("Creating minimum spanning tree...")
    sorted_weights_tuples = min_spanning_tree(knn_indices, knn_dists, n_neighbors, n_neighbors)
    
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
    
    return np.array(new_knn_indices), np.array(new_knn_dists)
