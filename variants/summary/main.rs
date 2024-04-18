use std::collections::HashMap;
use itertools::izip;
use hdf5;

mod helpers;
use helpers::*;

const QUAL_CUTOFF: u8 = 30; // 2 12 26 34

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        panic!("Usage: reads.h5");
    }
    let file = hdf5::File::open(&args[1]).expect("Could not open reads.h5 input file");
    
    // load the .h5 data
    let data = load_data(&file); // Vec<DataEntry> (read, rname_i, mapq, cb_i, ub_i)
    let order_map: HashMap<ReadNum, usize> = data.iter().map(|e| e.read).collect::<Vec<ReadNum>>().iter().enumerate().map(|(i, &num)| (num, i)).collect();
    let snv = load_snv(&file, &order_map); // Vec<SNVEntry> (read, pos, alt, qual)
    let matches = load_matches(&file, &order_map); // Vec<MatchEntry> (read, start, end)
    
    // Create data structures (umi) (cb_i) (ub_i) (reads) (avg_mapq)
    let mut data_umi: Vec<ReadNum> = Vec::new();
    let mut data_cb_i: Vec<WLi> = Vec::new();
    let mut data_ub_i: Vec<WLi> = Vec::new();
    let mut data_rname_i: Vec<RNAMEi> = Vec::new();
    let mut data_reads: Vec<ReadNum> = Vec::new();
    let mut data_mapq_avg: Vec<u8> = Vec::new();
    // Create snv table structure (snv_i) (RNAME) (POS) (ALT)
    let mut snv_table = SNVTable::new();
    // Create snv structures (umi) (snv_i) (hq) (lq) (total)
    let mut snv_umi: Vec<ReadNum> = Vec::new();
    let mut snv_snv_i: Vec<usize> = Vec::new();
    let mut snv_hq: Vec<ReadNum> = Vec::new();
    let mut snv_lq: Vec<ReadNum> = Vec::new();
    let mut snv_total: Vec<ReadNum> = Vec::new();
    // Create matches structures (umi) (start) (end)
    let mut match_umi: Vec<ReadNum> = Vec::new();
    let mut match_start: Vec<SeqPos> = Vec::new();
    let mut match_end: Vec<SeqPos> = Vec::new();
    // Create metadata structures
    let mut chimeric_reads: ReadNum = 0;
    
    println!("Processing reads...");

    // Loop through the reads
    let mut umi: ReadNum = 0;
    let mut curr_cb: WLi = data[0].cb_i;
    let mut curr_ub: WLi = data[0].ub_i;
    let mut snv_idx = 0; let mut matches_idx = 0;
    let mut umi_dict: HashMap<RNAMEi,UMIinfo> = HashMap::new(); // rname_i -> (reads, mapq, snv_i, start, end)
    for d in data { // d is (read, rname_i, mapq, cb_i, ub_i)
        if (d.cb_i != curr_cb) || (d.ub_i != curr_ub) {
          // remove chimerism, take dominant RNAME
          let total_reads: ReadNum = umi_dict.values().map(|record| record.reads).sum();
          let max_reads: ReadNum = umi_dict.values().map(|record| record.reads).max().expect("E");
          if umi_dict.values().filter(|record| record.reads == max_reads).count() > 1 { // if two RNAME tie for max read support, discard
              chimeric_reads += total_reads;
          } else { // get the RNAME with the most read support
              let (&rname_i, umi_info) = umi_dict.iter().max_by_key(|entry| entry.1.reads).expect("E");
              // write aggregate data
              data_umi.push(umi);
              data_cb_i.push(curr_cb); 
              data_ub_i.push(curr_ub);
              data_rname_i.push(rname_i);
              data_reads.push(umi_info.reads);
              data_mapq_avg.push(average_round_u8(&umi_info.mapq));
              // write aggregate snv
              for ((snv_i, pos), (hq, lq)) in umi_info.snv_i.iter() {
                  snv_umi.push(umi);
                  snv_snv_i.push(*snv_i);
                  snv_hq.push(*hq);
                  snv_lq.push(*lq);
                  snv_total.push(izip!(&umi_info.start, &umi_info.end).filter(|(&s, &e)| *pos >= s && *pos < e).count() as ReadNum);
              }
              // write aggregate matches
              let merged_intervals = merge_intervals(&umi_info.start, &umi_info.end);
              for interval in merged_intervals {
                match_umi.push(umi);
                match_start.push(interval.0);
                match_end.push(interval.1);
              }
              // write metadata
              chimeric_reads += total_reads - max_reads;
              umi += 1;
          }
          // reset umi struct
          umi_dict = HashMap::new();
          curr_cb = d.cb_i; curr_ub = d.ub_i;
        }
        let umi_info: &mut UMIinfo = umi_dict.entry(d.rname_i).or_insert(UMIinfo::new());
        umi_info.reads += 1;
        umi_info.mapq.push(d.mapq);

        // collect SNV data
        while snv_idx<snv.len() && snv[snv_idx].read==d.read { // snv is (read, pos, alt, qual)
            let snv_i: usize = snv_table.get(&(d.rname_i, snv[snv_idx].pos, snv[snv_idx].alt));
            let qual: u8 = snv[snv_idx].qual;
            umi_info.snv_i.entry((snv_i, snv[snv_idx].pos))
                .and_modify(|e| if qual > QUAL_CUTOFF { e.0 += 1 } else { e.1 += 1 })
                .or_insert(if qual > QUAL_CUTOFF { (1, 0) } else { (0, 1) });
            snv_idx += 1;
        }
        
        // collect matches data
        while matches_idx<matches.len() && matches[matches_idx].read==d.read { // matches is (read, start, end)
            umi_info.start.push(matches[matches_idx].start);
            umi_info.end.push(matches[matches_idx].end);
            matches_idx += 1;
        }
    }
    assert!(snv_idx == snv.len());
    assert!(matches_idx == matches.len());
    
    // let ins_read: Vec<u64> = file.dataset("ins/read").expect("not found").read_raw().expect("E");
    // let ins_pos: Vec<u64> = file.dataset("ins/pos").expect("not found").read_raw().expect("E"); assert!(ins_read.len() == ins_pos.len());
    // let ins_str_i: Vec<u64> = file.dataset("ins/str_i").expect("not found").read_raw().expect("E"); assert!(ins_read.len() == ins_str_i.len());
    
    // let del_read: Vec<u64> = file.dataset("del/read").expect("not found").read_raw().expect("E");
    // let del_pos: Vec<u64> = file.dataset("del/pos").expect("not found").read_raw().expect("E"); assert!(del_read.len() == del_pos.len());
    // let del_len: Vec<u64> = file.dataset("del/len").expect("not found").read_raw().expect("E"); assert!(del_read.len() == del_len.len());
    
    println!("Saving output...");
    
    // write umis data to an .h5 file
    let file = hdf5::File::create("umis.h5").expect("Could not create .h5 output file");
    // write general data
    let data_group = file.create_group("data").expect(".h5 error");
    assert_all_same(&[data_umi.len(), data_cb_i.len(), data_ub_i.len(), data_rname_i.len(), data_reads.len(), data_mapq_avg.len()]);
    data_umi.write_vector(&data_group, "umi");           // umi number (key)
    data_cb_i.write_vector(&data_group, "cb_i");         // cb whitelist index (0-indexed)
    data_ub_i.write_vector(&data_group, "ub_i");         // ub whitelist index (0-indexed)
    data_rname_i.write_vector(&data_group, "rname_i");   // rname whitelist index (0-indexed)
    data_reads.write_vector(&data_group, "reads");       // number of reads for the umi
    data_mapq_avg.write_vector(&data_group, "mapq_avg"); // integer average mapq for all reads of the umi
    // write the SNV table
    let tab_group = file.create_group("snv_table").expect(".h5 error");
    let (tab_snv_i, tab_rname_i, tab_pos, tab_alt) = snv_table.into_vectors();
    assert_all_same(&[tab_snv_i.len(), tab_rname_i.len(), tab_pos.len(), tab_alt.len()]);
    tab_snv_i.write_vector(&tab_group, "snv_i");     // SNV table index (0-indexed)
    tab_rname_i.write_vector(&tab_group, "rname_i"); // RNAME whitelist index (0-indexed)
    tab_pos.write_vector(&tab_group, "pos");         // POS of SNV
    tab_alt.write_vector(&tab_group, "alt");         // ALT base
    // write the SNV data
    let snv_group = file.create_group("snv").expect(".h5 error");
    assert_all_same(&[snv_umi.len(), snv_snv_i.len(), snv_hq.len(), snv_lq.len(), snv_total.len()]);
    snv_umi.write_vector(&snv_group, "umi");         // umi number (key)
    snv_snv_i.write_vector(&snv_group, "snv_i");     // SNV table index (0-indexed)
    snv_hq.write_vector(&snv_group, "hq");           // number of high-quality reads
    snv_lq.write_vector(&snv_group, "lq");           // number of low-quality reads
    snv_total.write_vector(&snv_group, "total");     // number of total reads
    // write the matches data
    let match_group = file.create_group("match").expect(".h5 error");
    assert_all_same(&[match_umi.len(), match_start.len(), match_end.len()]);
    match_umi.write_vector(&match_group, "umi");     // umi number (key)
    match_start.write_vector(&match_group, "start"); // 0-indexed position of range start (inclusive)
    match_end.write_vector(&match_group, "end");     // 0-indexed position of range end (exclusive)
    // write metadata
    let meta_group = file.create_group("metadata").expect(".h5 error");
    vec![chimeric_reads].write_vector(&meta_group, "reads_chimeric");
  
    println!("DONE");
}
