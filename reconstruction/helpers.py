from scipy.ndimage import gaussian_filter1d
from scipy.spatial import ConvexHull
from scipy.signal import find_peaks
from umi_tools import UMIClusterer
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import editdistance
import pandas as pd
import numpy as np

### matching bead barcode ###

def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km, mid_km, end_km

def barcode_matching(bc_pos_dict,spatial_bc_list,max_dist=1):
    bc_matching_dict = {}
    def get_sel_bc(bc):
        res = []
        if bc[:6] in start_km:
            res += start_km[bc[:6]]
        if bc[-6:] in end_km:
            res += end_km[bc[-6:]]
        if bc[4:10] in mid_km:
            res += mid_km[bc[4:10]]
        return set(res)
    exact_match = 0
    fuzzy_match = 0
    bc_ref_list = list(bc_pos_dict.keys())
    start_km,mid_km,end_km = build_6mer_dist(bc_ref_list)
    i=0
    for bc in spatial_bc_list:
        i+=1
        if i%500000==0:
            print(i)
        bc_old = bc
        #bc = "".join([bc[it] for it in [1,2,3,4,5,6,7,9,10,11,12,13] ])
        if bc in bc_pos_dict:
            exact_match += 1
            bc_matching_dict[bc_old] = bc
        else:
            sel_bc = get_sel_bc(bc)
            #sel_bc = bc_ref_list
            if len(sel_bc)>0:
                fz = [(it, editdistance.eval(it, bc)) for it in sel_bc]
                fz = [it for it in fz if it[1]<=max_dist]
                fz.sort(key=lambda x:x[1])
                if len(fz)==0:
                    continue
                ## if there are two barcodes with the same edit distance, choose the one with higher error rate in the last base
                if len(fz)>1 and fz[0][1]==fz[1][1]:
                    if editdistance.eval(fz[0][0][:-1], bc[-1])>editdistance.eval(fz[1][0][:-1], bc[-1]):  # higher error rate in the last base of the barcode
                        fuzzy_match += 1
                        bc_matching_dict[bc_old] = fz[1][0]
                    elif editdistance.eval(fz[0][0][:-1], bc[-1])<editdistance.eval(fz[1][0][:-1], bc[-1]):
                        fuzzy_match += 1
                        bc_matching_dict[bc_old] = fz[0][0]
                else:
                    fuzzy_match += 1
                    bc_matching_dict[bc_old] = fz[0][0]
    return bc_matching_dict, exact_match, fuzzy_match


### barcode collapsing ###

def umi_collapsing(cnt_dict, max_dist=1):
    """
    input: dict of barcode without matching
    output: list of barcode after collapsing
    """
    start_time = time.time()
    clusterer = UMIClusterer(cluster_method="directional")
    clustered_bc = clusterer(cnt_dict, threshold=max_dist)
    clustering_time = time.time()
    cluster_bc = [bc_group[0].decode('utf-8') for bc_group in clustered_bc]
    end_time = time.time()
    print("Clustering time: {}s".format(clustering_time-start_time))
    print("Dict creation time is: {}s".format(end_time-clustering_time))
    print("Total time is: {}s".format(end_time-start_time))
    return cluster_bc

def bc_collapsing(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1, min_reads_R2, alignment_stat):
    """ 
    input: dict of barcode without matching
    output: dict of barcode after filtering and collapsing
    """
    # filter for reads and collapse to whitelist
    R1_list = [s.encode('utf-8') for s in R1_bc_list]
    R1_dict = dict(Counter(R1_list))
    R1_dict_top = {k: v for k, v in R1_dict.items() if v > min_reads_R1}
    R1_whitelist = umi_collapsing(R1_dict_top)
    print("R1 total {}, after filter {}, whitelist {}".format(len(R1_dict),len(R1_dict_top),len(R1_whitelist)))
    print("read percentage: {}".format(np.sum(list(R1_dict_top.values()))/np.sum(list(R1_dict.values()))))
    R2_list = [s.encode('utf-8') for s in R2_bc_list]
    R2_dict = dict(Counter(R2_list))
    R2_dict_top = {k: v for k, v in R2_dict.items() if v > min_reads_R2}
    R2_whitelist = umi_collapsing(R2_dict_top)
    print("R2 total {}, after filter {}, whitelist {}".format(len(R2_dict),len(R2_dict_top),len(R2_whitelist)))
    print("read percentage: {}".format(np.sum(list(R2_dict_top.values()))/np.sum(list(R2_dict.values()))))

    # match to whitelist
    R1_bc_matching_dict,_,_ = barcode_matching(Counter(R1_whitelist), list(set(R1_bc_list)), max_dist=1)
    R2_bc_matching_dict,_,_ = barcode_matching(Counter(R2_whitelist), list(set(R2_bc_list)), max_dist=1)

    # generate dict with matched bc
    aln_dict_new = {}
    for bc_R1 in aln_dict:
        if bc_R1 in R1_bc_matching_dict:
            for R2 in range(len(aln_dict[bc_R1])):
                bc_R2 = aln_dict[bc_R1][R2][0]
                if bc_R2 in R2_bc_matching_dict:
                    alignment_stat["after_filter_reads"] += 1
                    aln_dict_new.setdefault(R1_bc_matching_dict[bc_R1],[]).append(
                        (R2_bc_matching_dict[bc_R2],aln_dict[bc_R1][R2][1],aln_dict[bc_R1][R2][2])) 
    print(len(aln_dict_new))
                    
    return aln_dict_new, alignment_stat
    
def bc_collecting(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1, min_reads_R2, alignment_stat):
    """ 
    input: dict of barcode without matching
    output: dict of barcode after filtering 
    """
    # filter for reads  to whitelist
    R1_dict = dict(Counter(R1_bc_list))
    R1_dict_top = {k: v for k, v in R1_dict.items() if v > min_reads_R1}
    R1_whitelist = set(R1_dict_top.keys())
    print("R1 total {}, after filter {}, whitelist {}".format(len(R1_dict),len(R1_dict_top),len(R1_whitelist)))
    print("read percentage: {}".format(np.sum(list(R1_dict_top.values()))/np.sum(list(R1_dict.values()))))
    R2_dict = dict(Counter(R2_bc_list))
    R2_dict_top = {k: v for k, v in R2_dict.items() if v > min_reads_R2}
    R2_whitelist = set(R2_dict_top.keys())
    print("R2 total {}, after filter {}, whitelist {}".format(len(R2_dict),len(R2_dict_top),len(R2_whitelist)))
    print("read percentage: {}".format(np.sum(list(R2_dict_top.values()))/np.sum(list(R2_dict.values()))))

    # generate dict with matched bc
    aln_dict_new = {}
    idx = 0
    for bc_R1 in aln_dict:
        idx += 1
        if bc_R1 in R1_whitelist:
            for R2 in range(len(aln_dict[bc_R1])):
                bc_R2 = aln_dict[bc_R1][R2][0]
                if bc_R2 in R2_whitelist:
                    alignment_stat["after_filter_reads"] += 1
                    aln_dict_new.setdefault(bc_R1,[]).append(
                        (bc_R2,aln_dict[bc_R1][R2][1],aln_dict[bc_R1][R2][2])) 
    print(len(aln_dict_new))
    return aln_dict_new, alignment_stat


### plots for fiducial_seq_blind_whitelist.py ###

def bc_rankplot(bc_list, position, qc_pdfs, max_expected_barcodes):
    bc_dict = Counter(bc_list).most_common()
    # sub = bc_dict[100:max_expected_barcodes].copy()
    sub = bc_dict
    x = np.histogram([np.log10(bc[1]) for bc in sub], bins=100)
    smooth = gaussian_filter1d(x[0], sigma=3)
    peak_idx,_ = find_peaks(-smooth)
    mean_hist = (x[1][1:][peak_idx]+x[1][:-1][peak_idx])/2
    print(f"{position} mean_hist: {mean_hist}")
    mean_hist = [m for m in mean_hist if m < 4]
    mean_hist = mean_hist[-1] # read cutoff

    bc_wl = [bc for bc in bc_dict if bc[1]>=10**mean_hist].copy()
    white_list_size=len(bc_wl)
    
    plt.figure(figsize=(4,3))
    log10_ranks=np.log10(np.arange(1,len(bc_dict)+1))
    log10_reads=[np.log10(bc[1]) for bc in bc_dict]
    plt.plot(log10_ranks,log10_reads)#,label='Rank Plot of Reads')
    plt.xlabel('Log10 Ranks')
    plt.ylabel('Log10 Reads')
    plt.title(f'{position} {white_list_size}')
    plt.plot([0, log10_ranks[-1]], [mean_hist, mean_hist], linewidth=1,label='log10 threshold',c='tab:green')
    log10_wl=np.log10(white_list_size)
    plt.plot([log10_wl, log10_wl], [np.min(log10_reads), np.max(log10_reads)], linewidth=1,label='log10 size',c='tab:orange')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    
    plt.figure(figsize=(4,3))
    plt.plot(x[1][0:-1], np.log10(x[0]+1), label='Raw Histogram')
    plt.plot(x[1][0:-1], np.log10(smooth+1), label='Gaussian Smoothed')
    plt.xlabel('Log10 UMI Counts')
    plt.ylabel('Log10 Bin Height')
    plt.title(f'{position}')
    plt.plot([mean_hist, mean_hist], [0, np.log10(np.max(x[0]+1))], linewidth=2, label='Whitelist Threshold')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    return 10**mean_hist


### plots for reconstruction_blind.py ###

# plot the distribution of umi each bead has
def blind_cnt_distribution(bead_all, bead_type, path):
    plt.figure(figsize=(8,6))
    sns.histplot(np.log10(bead_all['total_cnt']), bins=50)
    plt.xlabel('log10(total count)')
    plt.ylabel('number of '+bead_type)
    plt.title('blind '+bead_type+' total count distribution ({}), median={}'.format(len(bead_all), bead_all['total_cnt'].median()))
    plt.savefig(path, dpi=300)
    plt.close()

# plot the distribution of how many bc each bead covered
def blind_cover_bc_distribution(bead_cover_bc, bead_type, path):
    plt.figure(figsize=(8,6))
    sns.histplot(bead_cover_bc['cnt'], bins=50)
    plt.xlabel('bead covered')
    plt.ylabel('number of '+bead_type)
    plt.title(f'{bead_type} bead covered distribution ({len(bead_cover_bc)}), median={bead_cover_bc["cnt"].median()}')
    plt.savefig(path, dpi=300)
    plt.close()
    
# plot the shape 
def plot_shape(embedding, anchor, path):
    plt.figure(figsize=(6,6))
    plt.scatter(embedding[:,0], embedding[:,1], s=1)
    plt.title(f'{anchor} umap ({len(embedding)})')
    plt.savefig(path, dpi=300)
    plt.close()

# plot uniformness
def ot_random_ref(n):
    RADIUS = 1500
    r_vals = np.sqrt(np.random.uniform(size=n)) * RADIUS
    thet_vals = np.random.uniform(size=n) * 2 * np.pi
    a_random = np.column_stack((r_vals * np.cos(thet_vals), r_vals * np.sin(thet_vals)))
    a_random = a_random.astype(float)
    a_random = pd.DataFrame(a_random)
    a_random.columns = ['xcoord','ycoord']
    return a_random
def plot_uniformness(embedding, counts, path):
    a_random = ot_random_ref(counts.shape[0])
    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    im1 = sns.kdeplot(x=a_random.xcoord, y=a_random.ycoord, cmap="Blues", fill=True, levels=30, ax=axes[0])
    cbar1 = fig.colorbar(im1.collections[0], ax=axes[0])
    axes[0].set_title('random (max={:.7f})'.format(cbar1.vmax))
    im2 = sns.kdeplot(x=embedding[:,0],y=embedding[:,1], cmap="Blues", fill=True, levels=30, ax=axes[1])
    cbar2 = fig.colorbar(im2.collections[0], ax=axes[1])
    axes[1].set_title('recon (max={:.7f})'.format(cbar2.vmax))
    plt.savefig(path, dpi=300)
    plt.close()

# plot convex
def plot_convex(embedding, anchor, path):
    hull = ConvexHull(embedding)
    area = hull.volume
    perimeter = np.sum(np.linalg.norm(embedding[hull.vertices] - np.roll(embedding[hull.vertices], -1, axis=0), axis=1))
    plt.figure(figsize=(6,6))
    plt.scatter(embedding[:,0], embedding[:,1], s=1, alpha=0.8)
    for simplex in hull.simplices:
        plt.plot(embedding[simplex, 0], embedding[simplex, 1], 'r-', linewidth=1.5, alpha = 0.7)
    plt.title(anchor+' umap convex hull ({:.5f})'.format(perimeter**2/area / (4*np.pi)), fontsize=16, pad=20)
    plt.savefig(path, dpi=300)
    plt.close()

