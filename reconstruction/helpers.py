from scipy.ndimage import gaussian_filter1d
from scipy.spatial import ConvexHull
from scipy.signal import find_peaks
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import editdistance
import pandas as pd
import numpy as np

# matching bead barcode
def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km,mid_km,end_km
def barcode_matching(bc_pos_dict, spatial_bc_list, max_dist=1):
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
    return bc_matching_dict,exact_match,fuzzy_match


# Plotting functions (for fiducial_seq_blind_whitelist.py)

def bc_rankplot(bc_list, position, qc_pdfs, max_expected_barcodes=1000000):
    bc_dict = Counter(bc_list).most_common()
    sub = bc_dict[100:max_expected_barcodes].copy()
    x = np.histogram([np.log10(bc[1]) for bc in sub], 100)
    smooth = gaussian_filter1d(x[0], 3)
    peak_idx,_ = find_peaks(-smooth)
    if peak_idx is None:
        peak_idx = np.argmax(-smooth)
    mean_hist = (x[1][1:][peak_idx]+x[1][:-1][peak_idx])/2
    mean_hist = mean_hist[-1]

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
    plt.plot(x[1][:-1],x[0], label='Raw Histogram')
    plt.plot(x[1][:-1],smooth, label='Gaussian Smoothed')
    plt.xlabel('Log10 UMI Counts')
    plt.ylabel('Bin Height')
    plt.title(f'{position}')
    plt.plot([mean_hist, mean_hist], [0, np.max(x[0])], linewidth=2,label='Whitelist Threshold')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    return 10**mean_hist


# Plotting functions (for reconstruction_blind.py)

# plot the distribution of umi each bead has
def plot_blind_cnt_distribution(bead_all, bead_type, path):
    plt.figure(figsize=(8,6))
    sns.histplot(np.log10(bead_all['total_cnt']), bins=50)
    plt.xlabel('log10(total count)')
    plt.ylabel('number of '+bead_type)
    plt.title('blind '+bead_type+' total count distribution ({}), median={}'.format(len(bead_all), bead_all['total_cnt'].median()))
    plt.savefig(path, dpi=300)
    plt.close()

# plot the distribution of how many bc each bead covered
def plot_blind_cover_bc_distribution(bead_cover_bc, bead_type, path):
    plt.figure(figsize=(8,6))
    sns.histplot(bead_cover_bc['cnt'], bins=50)
    plt.xlabel('bead covered')
    plt.ylabel('number of '+bead_type)
    plt.title(bead_type+' bead covered distribution ({}), median={}'.format(len(bead_cover_bc), bead_cover_bc['cnt'].median()))
    plt.savefig(path, dpi=300)
    plt.close()

# plot the shape
def plot_umap(embedding, path):
    plt.figure(figsize=(6,6))
    plt.scatter(embedding[:,0], embedding[:,1], s=1)
    plt.title('umap ({})'.format(len(embedding)))
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
def plot_density(embedding, counts, path):
    a_random = ot_random_ref(len(counts))
    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    im1 = sns.kdeplot(x=a_random.xcoord, y=a_random.ycoord, cmap="Blues", fill=True, levels=30, ax=axes[0])
    cbar1 = fig.colorbar(im1.collections[0], ax=axes[0])
    axes[0].set_title('random (max={:.7f})'.format(cbar1.vmax))
    im2 = sns.kdeplot(x=embedding[:,0],y=embedding[:,1], cmap="Blues", fill=True, levels=30, ax=axes[1])
    cbar2 = fig.colorbar(im2.collections[0], ax=axes[1])
    axes[1].set_title('recon (max={:.7f})'.format(cbar2.vmax))
    plt.savefig(path, dpi=300)
    plt.close()

# plot covex
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
