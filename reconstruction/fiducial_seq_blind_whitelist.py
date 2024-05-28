"""
analysis of fiducial diffusion sequencing result without bead barcode matching
input fastq file
output collasped barcode information
"""
import os
import time
import gzip
import argparse
import mappy as mp
import editdistance
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
from helpers import *

def barcode_extract(fq1_file, fq2_file, r2type):
    """
    input: fastq file
    output: dict of barcode without matching
    """
    aln_dict = {}
    R1_bc_list = []
    R2_bc_list = []
    alignment_stat = Counter()
    for fq1, fq2 in zip(mp.fastx_read(fq1_file, read_comment=False), mp.fastx_read(fq2_file, read_comment=False)):
        alignment_stat["total_reads"] += 1

        if alignment_stat["total_reads"] % 1000000 == 0:
            print(alignment_stat["total_reads"])
    
        # check read length
        if len(fq1[1]) < 46 or len(fq2[1]) < 46:  
            alignment_stat["Read_too_short"] += 1
            continue
        R1_bc = fq1[1][0:8] + fq1[1][26:33]
        R1_bumi = fq1[1][33:42] # if V5: 32:40, if V8 or V10: 33:41
        R1_UP = fq1[1][8:26]

        # different Read 2 bead type
        if r2type == 'V9':
            R2_bc = fq2[1][0:8] + fq2[1][26:33]
            R2_bumi = fq2[1][33:42]
            R2_UP = fq2[1][8:26]
            if editdistance.eval(R1_UP,'TCTTCAGCGTTCCCGAGA')>3 or editdistance.eval(R2_UP,'TCTTCAGCGTTCCCGAGA')>3:
                alignment_stat["UP_not_matched"] += 1
                continue 
        elif r2type == 'V15':
            R2_bc = fq2[1][0:15]
            R2_bumi = fq2[1][25:34]
            R2_UP = fq2[1][15:25]
            if editdistance.eval(R1_UP,'TCTTCAGCGTTCCCGAGA')>3 or editdistance.eval(R2_UP,'CTGTTTCCTG')>2:
                alignment_stat["UP_not_matched"] += 1
                continue
        else:
            assert False, f"unknown read2type ({r2type})"
          
        alignment_stat['effective_read'] += 1
        aln_dict.setdefault(R1_bc, []).append((R2_bc, R1_bumi, R2_bumi)) 
        R1_bc_list.append(R1_bc)
        R2_bc_list.append(R2_bc)
    return aln_dict, alignment_stat, R1_bc_list, R2_bc_list

def write_blind(aln_dict_new, alignment_stat, out_dir):
    # collapse for reads
    for bc in aln_dict_new:
        tmp = Counter(aln_dict_new[bc])
        aln_dict_new[bc] = tmp

    # write result to csv
    raw_f = gzip.open(os.path.join(out_dir,"blind_raw_reads_filtered.csv.gz"),"wb")
    raw_f.write(b'R1_bc,R2_bc,R1_bumi,R2_bumi,reads\n')
    for bc_R1 in aln_dict_new:
        raw_f.write(bytes('\n'.join(['{},{},{},{},{}'.format(
            bc_R1, it[0], it[1], it[2], aln_dict_new[bc_R1][it]) for it in aln_dict_new[bc_R1]])+'\n',"UTF-8"))
    raw_f.close()
    print("Write matched data to {}".format("blind_raw_reads_filtered.csv.gz"))

    with open(os.path.join(out_dir,"blind_statistics_filtered.csv"),"w") as f:
        f.write("alignment_status,counts\n")
        for aln_stat in alignment_stat:
            f.write(f"{aln_stat},{alignment_stat[aln_stat]}\n")

def get_args():
    parser = argparse.ArgumentParser(description='Process recon seq data.')
    parser.add_argument("-i", "--in_dir",
        help="input data folder",
        type=str,
        required=True,
    )
    parser.add_argument("-o", "--out_dir",
        help="output data folder",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r", "--read2type",
        help="input bead type of read2 (V9 or V15)",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    in_dir = args.in_dir
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    read2type = args.read2type

    print("loading files")
    fq_files = [f for f in os.listdir(in_dir) if not f.startswith('.')]
    R1s = [it for it in fq_files if "_R1_" in it]
    R2s = [it for it in fq_files if "_R2_" in it]
    print(f"R1s: {R1s}")
    print(f"R2s: {R2s}")
    assert len(R1s) == len(R2s) == 1
    fq1_file = os.path.join(in_dir, R1s[0])
    fq2_file = os.path.join(in_dir, R2s[0])
    
    # barode blind
    print("extracting barcodes")
    aln_dict, stat, R1_bc_list, R2_bc_list = barcode_extract(fq1_file, fq2_file, read2type)
    
    print("creating barcode rank plots")
    qc_pdf_file = os.path.join(out_dir, 'QC.pdf')
    qc_pdfs = PdfPages(qc_pdf_file)
    R1_threshold = bc_rankplot(R1_bc_list, 'R1', qc_pdfs, max_expected_barcodes=1000000)
    R2_threshold = bc_rankplot(R2_bc_list, 'R2', qc_pdfs, max_expected_barcodes=1000000)
    qc_pdfs.close()

    # R1_threshold = 50
    # R2_threshold = 50
    print(f"R1_threshold: {R1_threshold}")
    print(f"R2_threshold: {R2_threshold}")

    print("performing barcode collapsing")
    # if do barcode collapsing
    # aln_dict_new, stat_new = bc_collapsing(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1=R1_threshold, min_reads_R2=R2_threshold, alignment_stat = stat)
    # if no barcode collapsing
    aln_dict_new, stat_new = bc_collecting(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1=R1_threshold, min_reads_R2=R2_threshold, alignment_stat=stat)
    
    print("writing results")
    write_blind(aln_dict_new, stat_new, out_dir)
