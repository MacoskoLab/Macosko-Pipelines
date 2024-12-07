use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use hdf5;

use arrow::record_batch::RecordBatch;
use arrow::array::{UInt64Array, UInt32Array, UInt16Array, UInt8Array};
use arrow::array::ArrayRef;
use std::sync::Arc;
use arrow::datatypes::Schema;
use parquet::arrow::ArrowWriter;

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
    
    // PROCESS READS //
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

    // Loop through the reads
    let mut umi: ReadNum = 0;
    let mut snv_idx: usize = 0;
    let mut matches_idx: usize = 0;
    let mut umi_dict: HashMap<(RNAMEi, u8), UMIinfo> = HashMap::new(); // (rname_i, strand) -> (reads, mapq, snv, matches)
    for (data_idx, d) in data.iter().enumerate() { // d is (read, flag, rname_i, mapq, cb_i, ub_i)
    
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
        
        // check if there are more reads to aggregate
        if data_idx+1 < data.len() && d.cb_i == data[data_idx+1].cb_i && d.ub_i == data[data_idx+1].ub_i {
            continue
        }
        
        // if two (RNAME, strand) tie for max read support, discard the whole umi
        let total_reads: ReadNum = umi_dict.values().map(|record| record.reads).sum();
        let max_reads: ReadNum = umi_dict.values().map(|record| record.reads).max().unwrap();
        if umi_dict.values().filter(|ui| ui.reads == max_reads).count() > 1 { 
            chimeric_reads += total_reads;
            chimeric_umis += 1;
            umi_dict = HashMap::new();
            continue
        }
        
        // save the (RNAME, strand) with the most read support
        let (&(rname_i, strand), umi_info) = umi_dict.iter().max_by_key(|entry| entry.1.reads).unwrap();
        umi += 1;
        
        // write aggregate data
        data_umi.push(umi);
        data_cb_i.push(d.cb_i); 
        data_ub_i.push(d.ub_i);
        data_rname_i.push(rname_i);
        data_strand.push(strand);
        data_reads.push(umi_info.reads);
        data_mapq_avg.push(umi_info.mapq_avg());
        
        // write aggregate snv
        for (&(pos, alt), &(hq, lq)) in umi_info.snv.iter() {
            snv_umi.push(umi);
            snv_pos.push(pos);
            snv_alt.push(alt);
            snv_hq.push(hq);
            snv_lq.push(lq);
            snv_total.push(umi_info.matches.iter()
                                           .filter(|(s, e)| pos >= *s && pos < *e)
                                           .count()
                                           .try_into()
                                           .expect("snv_total data type too small"));
        }
        
        // write aggregate matches
        for interval in umi_info.merge_intervals() {
          match_umi.push(umi);
          match_start.push(interval.0);
          match_end.push(interval.1);
        }
        
        // update metadata
        chimeric_reads += total_reads - umi_info.reads;
        
        // reset umi struct
        umi_dict = HashMap::new();
    }
    assert!(snv_idx == snv.len());
    assert!(matches_idx == matches.len());
    
    // SAVE RESULTS //
    eprintln!("Saving output...");
    
    // Create .parquet files
    let out_dir = Path::new(&args[1]).parent().unwrap();
    let data_file = File::create(&out_dir.join("data.parquet")).expect("Could not create data.parquet output file");
    let snv_file = File::create(&out_dir.join("snv.parquet")).expect("Could not create snv.parquet output file");
    let match_file = File::create(&out_dir.join("match.parquet")).expect("Could not create match.parquet output file");
    
    // Pack metadata
    let metadata: HashMap<String, String> = HashMap::from([
        ("chimeric_reads".to_string(), chimeric_reads.to_string()),
        ("chimeric_umis".to_string(), chimeric_umis.to_string())
    ]);
    
    // Write general data
    let data_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(data_umi)) as ArrayRef),          // umi number (key)
        ("cb_i", Arc::new(UInt32Array::from(data_cb_i)) as ArrayRef),        // cb whitelist index (0-indexed)
        ("ub_i", Arc::new(UInt32Array::from(data_ub_i)) as ArrayRef),        // ub whitelist index (0-indexed)
        ("rname_i", Arc::new(UInt16Array::from(data_rname_i)) as ArrayRef),  // rname whitelist index (0-indexed)
        ("strand", Arc::new(UInt8Array::from(data_strand)) as ArrayRef),     // 0 normal, 1 means reverse complemented
        ("reads", Arc::new(UInt64Array::from(data_reads)) as ArrayRef),      // number of reads for the umi
        ("mapq_avg", Arc::new(UInt8Array::from(data_mapq_avg)) as ArrayRef), // integer average mapq for all reads of the umi
    ]).unwrap();
    let data_schema = Schema::new_with_metadata(data_batch.schema().fields().clone(), metadata);
    let mut writer = ArrowWriter::try_new(data_file, Arc::new(data_schema), None).unwrap();
    writer.write(&data_batch).unwrap();
    writer.close().unwrap();
    
    // Write the snv data
    let snv_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(snv_umi)) as ArrayRef),     // umi number (key)
        ("pos", Arc::new(UInt32Array::from(snv_pos)) as ArrayRef),     // POS of SNV
        ("alt", Arc::new(UInt8Array::from(snv_alt)) as ArrayRef),      // ALT base
        ("hq", Arc::new(UInt64Array::from(snv_hq)) as ArrayRef),       // number of high-quality reads
        ("lq", Arc::new(UInt64Array::from(snv_lq)) as ArrayRef),       // number of low-quality reads
        ("total", Arc::new(UInt64Array::from(snv_total)) as ArrayRef), // number of total reads observing its POS
    ]).unwrap();
    let mut writer = ArrowWriter::try_new(snv_file, snv_batch.schema(), None).unwrap();
    writer.write(&snv_batch).unwrap();
    writer.close().unwrap();
    
    // Write the matches data
    let match_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(match_umi)) as ArrayRef),     // umi number (key)
        ("start", Arc::new(UInt32Array::from(match_start)) as ArrayRef), // 0-indexed position of range start (inclusive)
        ("end", Arc::new(UInt32Array::from(match_end)) as ArrayRef),     // 0-indexed position of range end (exclusive)
    ]).unwrap();
    let mut writer = ArrowWriter::try_new(match_file, match_batch.schema(), None).unwrap();
    writer.write(&match_batch).unwrap();
    writer.close().unwrap();
    
    eprintln!("Done!");
}
