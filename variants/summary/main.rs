use std::collections::HashMap;
use std::path::Path;
use hdf5;

mod helpers;
use helpers::*;

const QUAL_CUTOFF: u8 = 30; // 2 12 26 34

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        panic!("Usage: summary reads.h5");
    }
    
    eprintln!("Loading reads...");
    
    let file = hdf5::File::open(&args[1]).expect("Could not open reads.h5 input file");
    
    // Load the .h5 data
    let mut data = load_data(&file);       // Vec<DataEntry>  (read, flag, rname_i, mapq, cb_i, ub_i)
    let mut snv = load_snv(&file);         // Vec<SNVEntry>   (read, pos, alt, qual)
    let mut matches = load_matches(&file); // Vec<MatchEntry> (read, start, end)
    
    eprintln!("Sorting reads...");
    
    // Sort the .h5 data
    data.sort_by(|d1, d2| { (d1.cb_i).cmp(&d2.cb_i).then((d1.ub_i).cmp(&d2.ub_i)) }); // Sort by cb, umi
    let order_map: HashMap<ReadNum, usize> = data.iter().enumerate().map(|(i, e)| (e.read, i)).collect();
    snv.sort_by_key(|e| order_map.get(&e.read).unwrap());
    matches.sort_by_key(|e| order_map.get(&e.read).unwrap());
    
    eprintln!("Processing reads...");
    
    // Create data structures (umi) (cb_i) (ub_i) (rname_i) (strand) (reads) (mapq_avg)
    let mut data_umi: Vec<ReadNum> = Vec::new();
    let mut data_cb_i: Vec<WLi> = Vec::new();
    let mut data_ub_i: Vec<WLi> = Vec::new();
    let mut data_rname_i: Vec<RNAMEi> = Vec::new();
    let mut data_strand: Vec<u8> = Vec::new();
    let mut data_reads: Vec<ReadNum> = Vec::new();
    let mut data_mapq_avg: Vec<u8> = Vec::new();
    // Create snv structures (umi) (pos) (alt) (hq) (lq) (total)
    let mut snv_umi: Vec<ReadNum> = Vec::new();
    let mut snv_pos: Vec<SeqPos> = Vec::new();
    let mut snv_alt: Vec<u8> = Vec::new();
    let mut snv_hq: Vec<ReadNum> = Vec::new();
    let mut snv_lq: Vec<ReadNum> = Vec::new();
    let mut snv_total: Vec<ReadNum> = Vec::new();
    // Create matches structures (umi) (start) (end)
    let mut match_umi: Vec<ReadNum> = Vec::new();
    let mut match_start: Vec<SeqPos> = Vec::new();
    let mut match_end: Vec<SeqPos> = Vec::new();
    // Create metadata structures
    let mut chimeric_reads: ReadNum = 0;
    let mut chimeric_umis: ReadNum = 0;

    // Initial state
    let mut curr_cb: WLi = data.get(0).expect("ERROR: empty").cb_i;
    let mut curr_ub: WLi = data.get(0).expect("ERROR: empty").ub_i;
    let mut umi_dict: HashMap<(RNAMEi, u8), UMIinfo> = HashMap::new(); // (rname_i, strand) -> (reads, mapq, snv, matches)
    // Loop through the reads
    let mut umi: ReadNum = 0;
    let mut snv_idx: usize = 0;
    let mut matches_idx: usize = 0;
    for d in data { // d is (read, flag, rname_i, mapq, cb_i, ub_i)
    
        if (d.cb_i != curr_cb) || (d.ub_i != curr_ub) {
          // remove chimerism, take dominant RNAME
          let total_reads: ReadNum = umi_dict.values().map(|record| record.reads).sum();
          let max_reads: ReadNum = umi_dict.values().map(|record| record.reads).max().unwrap();
          if umi_dict.values().filter(|record| record.reads == max_reads).count() > 1 { // if two RNAME tie for max read support, discard the whole umi
              chimeric_reads += total_reads;
              chimeric_umis += 1;
          } else { // get the RNAME with the most read support
              let (&(rname_i, strand), umi_info) = umi_dict.iter().max_by_key(|entry| entry.1.reads).unwrap();
              // write aggregate data
              data_umi.push(umi);
              data_cb_i.push(curr_cb); 
              data_ub_i.push(curr_ub);
              data_rname_i.push(rname_i);
              data_strand.push(strand);
              data_reads.push(umi_info.reads);
              data_mapq_avg.push(umi_info.mapq_avg());
              // write aggregate snv
              for ((pos, alt), (hq, lq)) in umi_info.snv.iter() {
                  snv_umi.push(umi);
                  snv_pos.push(*pos);
                  snv_alt.push(*alt);
                  snv_hq.push(*hq);
                  snv_lq.push(*lq);
                  snv_total.push(umi_info.matches.iter().filter(|&(s, e)| *pos >= *s && *pos < *e).count() as ReadNum);
              }
              // write aggregate matches
              for interval in umi_info.merge_intervals() {
                match_umi.push(umi);
                match_start.push(interval.0);
                match_end.push(interval.1);
              }
              // write metadata
              chimeric_reads += total_reads - umi_info.reads;
              umi += 1;
          }
          // reset umi struct
          umi_dict = HashMap::new();
          curr_cb = d.cb_i;
          curr_ub = d.ub_i;
        }
        
        // get the UMIinfo for the current (rname_i, strand)
        let umi_info: &mut UMIinfo = umi_dict.entry((d.rname_i, flag2strand(d.flag))).or_insert(UMIinfo::new()); 
        
        // collect general data
        umi_info.reads += 1;
        umi_info.mapq.push(d.mapq);

        // collect SNV data
        while snv_idx < snv.len() && snv[snv_idx].read == d.read { // snv is (read, pos, alt, qual)
            let pos = snv[snv_idx].pos;
            let alt = snv[snv_idx].alt;
            let qual = snv[snv_idx].qual;
            umi_info.snv.entry((pos, alt))
                        .and_modify(|e| if qual > QUAL_CUTOFF { e.0 += 1 } else { e.1 += 1 })
                        .or_insert(if qual > QUAL_CUTOFF { (1, 0) } else { (0, 1) });
            snv_idx += 1;
        }
        
        // collect matches data
        while matches_idx < matches.len() && matches[matches_idx].read == d.read { // matches is (read, start, end)
            let start = matches[matches_idx].start;
            let end = matches[matches_idx].end;
            umi_info.matches.push((start, end));
            matches_idx += 1;
        }
    }
    assert!(snv_idx == snv.len());
    assert!(matches_idx == matches.len());
    
    eprintln!("Saving output...");
    
    // Create umis.h5 file
    let out_dir = Path::new(&args[1]).parent().unwrap();
    let file = hdf5::File::create(out_dir.join("umis.h5")).expect("Could not create .h5 output file");
    // Write general data
    let data_group = file.create_group("data").expect(".h5 error");
    assert_all_same(&[data_umi.len(), data_cb_i.len(), data_ub_i.len(), data_rname_i.len(), data_strand.len(), data_reads.len(), data_mapq_avg.len()]);
    data_umi.write_vector(&data_group, "umi");           // umi number (key)
    data_cb_i.write_vector(&data_group, "cb_i");         // cb whitelist index (0-indexed)
    data_ub_i.write_vector(&data_group, "ub_i");         // ub whitelist index (0-indexed)
    data_rname_i.write_vector(&data_group, "rname_i");   // rname whitelist index (0-indexed)
    data_strand.write_vector(&data_group, "strand");     // 0 normal, 1 means reverse complemented
    data_reads.write_vector(&data_group, "reads");       // number of reads for the umi
    data_mapq_avg.write_vector(&data_group, "mapq_avg"); // integer average mapq for all reads of the umi
    // Write the SNV data
    let snv_group = file.create_group("snv").expect(".h5 error");
    assert_all_same(&[snv_umi.len(), snv_pos.len(), snv_alt.len(), snv_hq.len(), snv_lq.len(), snv_total.len()]);
    snv_umi.write_vector(&snv_group, "umi");         // umi number (key)
    snv_pos.write_vector(&snv_group, "pos");         // POS of SNV
    snv_alt.write_vector(&snv_group, "alt");         // ALT base
    snv_hq.write_vector(&snv_group, "hq");           // number of high-quality reads
    snv_lq.write_vector(&snv_group, "lq");           // number of low-quality reads
    snv_total.write_vector(&snv_group, "total");     // number of total reads observing its POS
    // Write the matches data
    let match_group = file.create_group("match").expect(".h5 error");
    assert_all_same(&[match_umi.len(), match_start.len(), match_end.len()]);
    match_umi.write_vector(&match_group, "umi");     // umi number (key)
    match_start.write_vector(&match_group, "start"); // 0-indexed position of range start (inclusive)
    match_end.write_vector(&match_group, "end");     // 0-indexed position of range end (exclusive)
    // Write metadata
    let meta_group = file.create_group("metadata").expect(".h5 error");
    vec![chimeric_reads].write_vector(&meta_group, "reads_chimeric");
    vec![chimeric_umis].write_vector(&meta_group, "umis_chimeric");
    
    eprintln!("Done!");
}
