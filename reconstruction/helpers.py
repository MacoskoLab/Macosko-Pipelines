import os
from collections import Counter
import editdistance

### matching bead barcode ###
def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km,mid_km,end_km
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
    return bc_matching_dict,exact_match,fuzzy_match

### plots for reconstruction_blind.py ###

# plot the distribution of umi each bead has
def blind_cnt_distribution(bead_all, bead_type, path=os.path.join(out_dir,bead_type+'_blind_cnt_distribution.png')):
    plt.figure(figsize=(8,6))
    sns.histplot(np.log10(bead_all['total_cnt']), bins=50)
    plt.xlabel('log10(total count)')
    plt.ylabel('number of '+bead_type)
    plt.title('blind '+bead_type+' total count distribution ({}), median={}'.format(len(bead_all), bead_all['total_cnt'].median()))
    plt.savefig(path, dpi=300)
    plt.close()

# plot the distribution of how many bc each bead covered
def blind_cover_bc_distribution(bead_cover_bc, bead_type, path=os.path.join(out_dir,bead_type+'_blind_cover_bc_distribution.png')):
    plt.figure(figsize=(8,6))
    sns.histplot(bead_cover_bc['cnt'], bins=50)
    plt.xlabel('bead covered')
    plt.ylabel('number of '+bead_type)
    plt.title(bead_type+' bead covered distribution ({}), median={}'.format(len(bead_cover_bc), bead_cover_bc['cnt'].median()))
    plt.savefig(path, dpi=300)
    plt.close()
