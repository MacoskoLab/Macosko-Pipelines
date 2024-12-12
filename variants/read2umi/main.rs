use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use hdf5;

use arrow::record_batch::RecordBatch;
use arrow::array::{UInt64Array, UInt32Array, UInt16Array, UInt8Array};
use arrow::array::{UInt32DictionaryArray, UInt16DictionaryArray};
use arrow::array::StringArray;
use arrow::array::ArrayRef;
use std::sync::Arc;
use arrow::datatypes::Schema;
use parquet::arrow::ArrowWriter;

mod helpers;
use helpers::*;

const QUAL_CUTOFF: u8 = 30; // 2 12 26 34
const BUFF_CUTOFF: SeqPos = 15;

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
    let mut ins = load_insertions(&file);  // Vec<InsEntry> (read, pos, str_i)
    let mut del = load_deletions(&file);   // Vec<DelEntry> (read, pos, len)
    
    eprintln!("Sorting reads...");
    
    // Sort the .h5 data
    data.sort_by(|d1, d2| { (d1.cb_i).cmp(&d2.cb_i).then((d1.ub_i).cmp(&d2.ub_i)) }); // Sort by cb, umi
    let order_map: HashMap<ReadNum, usize> = data.iter().enumerate().map(|(i, e)| (e.read, i)).collect();
    snv.sort_by_key(|e| order_map.get(&e.read).unwrap());
    matches.sort_by_key(|e| order_map.get(&e.read).unwrap());
    ins.sort_by_key(|e| order_map.get(&e.read).unwrap());
    del.sort_by_key(|e| order_map.get(&e.read).unwrap());
    
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
    // Create snv structures (umi) (pos) (alt) (hq) (lq) (covers)
    let mut snv_umi: Vec<ReadNum> = Vec::new();
    let mut snv_pos: Vec<SeqPos> = Vec::new();
    let mut snv_alt: Vec<u8> = Vec::new();
    let mut snv_hq: Vec<ReadNum> = Vec::new();
    let mut snv_lq: Vec<ReadNum> = Vec::new();
    let mut snv_covers: Vec<ReadNum> = Vec::new();
    // Create matches structures (umi) (start) (end)
    let mut match_umi: Vec<ReadNum> = Vec::new();
    let mut match_start: Vec<SeqPos> = Vec::new();
    let mut match_end: Vec<SeqPos> = Vec::new();
    // Create insertion structures
    let mut ins_umi: Vec<ReadNum> = Vec::new();
    let mut ins_pos: Vec<SeqPos> = Vec::new();
    let mut ins_str_i: Vec<WLi> = Vec::new();
    let mut ins_obs: Vec<ReadNum> = Vec::new();
    let mut ins_hq: Vec<ReadNum> = Vec::new();
    let mut ins_lq: Vec<ReadNum> = Vec::new();
    // Create deletion structures
    let mut del_umi: Vec<ReadNum> = Vec::new();
    let mut del_pos: Vec<SeqPos> = Vec::new();
    let mut del_len: Vec<SeqPos> = Vec::new();
    let mut del_obs: Vec<ReadNum> = Vec::new();
    let mut del_hq: Vec<ReadNum> = Vec::new();
    let mut del_lq: Vec<ReadNum> = Vec::new();
    // Create metadata structures
    let mut chimeric_reads: ReadNum = 0;
    let mut chimeric_umis: ReadNum = 0;

    // Loop through the reads
    let mut umi: ReadNum = 0;
    let mut snv_idx: usize = 0;
    let mut matches_idx: usize = 0;
    let mut ins_idx: usize = 0;
    let mut del_idx: usize = 0;
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
        
        // collect insertion data
        while ins_idx < ins.len() && ins[ins_idx].read == d.read { // ins is (read, pos, str_i)
            let pos = ins[ins_idx].pos;
            let str_i = ins[ins_idx].str_i;
            *umi_info.ins.entry((pos, str_i)).or_insert(0) += 1;
            ins_idx += 1;
        }
        
        // collect deletion data
        while del_idx < del.len() && del[del_idx].read == d.read { // del is (read, pos, len)
            let pos = del[del_idx].pos;
            let len = del[del_idx].len;
            *umi_info.del.entry((pos, len)).or_insert(0) += 1;
            del_idx += 1;
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
            snv_covers.push(umi_info.matches.iter()
                                            .filter(|(s, e)| pos >= *s && pos < *e)
                                            .count().try_into().expect("snv_covers data type too small"));
        }
        
        // write aggregate matches
        for interval in umi_info.merge_intervals() {
            match_umi.push(umi);
            match_start.push(interval.0);
            match_end.push(interval.1);
        }
        
        // write aggregate insertions
        for (&(pos, str_i), &obs) in umi_info.ins.iter() {
            ins_umi.push(umi);
            ins_pos.push(pos);
            ins_str_i.push(str_i);
            ins_obs.push(obs);
            let covers_all: ReadNum = umi_info.matches.iter()
                                              .filter(|(s, e)| *s < pos && *e > pos)
                                              .count().try_into().unwrap();
            let covers_hq: ReadNum = umi_info.matches.iter()
                                             .filter(|(s, e)| *s + BUFF_CUTOFF < pos && *e - BUFF_CUTOFF > pos)
                                             .count().try_into().unwrap();
            ins_hq.push(covers_hq);
            ins_lq.push(covers_all - covers_hq);
        }
        
        // write aggregate deletions
        for (&(pos, len), &obs) in umi_info.del.iter() {
            del_umi.push(umi);
            del_pos.push(pos);
            del_len.push(len);
            del_obs.push(obs);
            let covers_all: ReadNum = umi_info.matches.iter()
                                              .filter(|(s, e)| *s < pos && *e > pos)
                                              .count().try_into().unwrap();
            let covers_hq: ReadNum = umi_info.matches.iter()
                                             .filter(|(s, e)| *s + BUFF_CUTOFF < pos && *e - BUFF_CUTOFF > pos)
                                             .count().try_into().unwrap();
            del_hq.push(covers_hq);
            del_lq.push(covers_all - covers_hq);
        }
        
        // update metadata
        chimeric_reads += total_reads - umi_info.reads;
        
        // reset umi struct
        umi_dict = HashMap::new();
    }
    assert!(snv_idx == snv.len());
    assert!(matches_idx == matches.len());
    assert!(ins_idx == ins.len());
    assert!(del_idx == del.len());
    
    std::mem::drop(data);
    std::mem::drop(snv);
    std::mem::drop(matches);
    std::mem::drop(ins);
    std::mem::drop(del);
    
    // SAVE RESULTS //
    eprintln!("Saving output...");
    
    // Turn index+whitelist into arrow dictionary
    use hdf5::types::VarLenAscii;
    let cb: Vec<VarLenAscii> = file.dataset("whitelists/cb").expect("not found").read_raw().unwrap();
    let ub: Vec<VarLenAscii> = file.dataset("whitelists/ub").expect("not found").read_raw().unwrap();
    let ins: Vec<VarLenAscii> = file.dataset("whitelists/ins").expect("not found").read_raw().unwrap();
    // "whitelists/sc"
    let rname: Vec<VarLenAscii> = file.dataset("whitelists/rname").expect("not found").read_raw().unwrap();
    
    let data_cb = UInt32DictionaryArray::new(UInt32Array::from(data_cb_i), Arc::new(StringArray::from_iter_values(cb)) as ArrayRef);
    let data_ub = UInt32DictionaryArray::new(UInt32Array::from(data_ub_i), Arc::new(StringArray::from_iter_values(ub)) as ArrayRef);
    let data_rname = UInt16DictionaryArray::new(UInt16Array::from(data_rname_i), Arc::new(StringArray::from_iter_values(rname)) as ArrayRef);
    let ins_str = UInt32DictionaryArray::new(UInt32Array::from(ins_str_i), Arc::new(StringArray::from_iter_values(ins)) as ArrayRef);
    
    // Pack metadata
    let metadata: HashMap<String, String> = HashMap::from([
        ("chimeric_reads".to_string(), chimeric_reads.to_string()),
        ("chimeric_umis".to_string(), chimeric_umis.to_string())
    ]);
    
    // Create .parquet files
    let out_dir = Path::new(&args[1]).parent().unwrap();
    let data_file = File::create(&out_dir.join("data.parquet")).expect("Could not create data.parquet output file");
    let snv_file = File::create(&out_dir.join("snv.parquet")).expect("Could not create snv.parquet output file");
    let match_file = File::create(&out_dir.join("match.parquet")).expect("Could not create match.parquet output file");
    let ins_file = File::create(&out_dir.join("ins.parquet")).expect("Could not create ins.parquet output file");
    let del_file = File::create(&out_dir.join("del.parquet")).expect("Could not create del.parquet output file");
    
    // Write general data
    let data_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(data_umi)) as ArrayRef),          // umi number (key)
        ("cb", Arc::new(data_cb) as ArrayRef),                               // cell barcode
        ("ub", Arc::new(data_ub) as ArrayRef),                               // UMI barcode
        ("rname", Arc::new(data_rname) as ArrayRef),                         // RNAME
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
        ("umi", Arc::new(UInt64Array::from(snv_umi)) as ArrayRef),       // umi number (key)
        ("pos", Arc::new(UInt32Array::from(snv_pos)) as ArrayRef),       // 0-indexed position of SNV
        ("alt", Arc::new(UInt8Array::from(snv_alt)) as ArrayRef),        // ALT base
        ("hq", Arc::new(UInt64Array::from(snv_hq)) as ArrayRef),         // number of reads with high-quality alt base
        ("lq", Arc::new(UInt64Array::from(snv_lq)) as ArrayRef),         // number of reads with low-quality alt base
        ("covers", Arc::new(UInt64Array::from(snv_covers)) as ArrayRef), // number of reads observing its pos
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
    
    // Write the insertion data
    let ins_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(ins_umi)) as ArrayRef), // umi number (key)
        ("pos", Arc::new(UInt32Array::from(ins_pos)) as ArrayRef), // 0-indexed position of reference base after insertion
        ("str", Arc::new(ins_str) as ArrayRef),                    // inserted string
        ("obs", Arc::new(UInt64Array::from(ins_obs)) as ArrayRef), // number of reads with this insertion
        ("hq", Arc::new(UInt64Array::from(ins_hq)) as ArrayRef),   // number of reads that did NOT observe the insertion (high confidence)
        ("lq", Arc::new(UInt64Array::from(ins_lq)) as ArrayRef),   // number of reads that did NOT observe the insertion (low confidence)
    ]).unwrap();
    let mut writer = ArrowWriter::try_new(ins_file, ins_batch.schema(), None).unwrap();
    writer.write(&ins_batch).unwrap();
    writer.close().unwrap();
    
    // Write the deletion data
    let del_batch = RecordBatch::try_from_iter(vec![
        ("umi", Arc::new(UInt64Array::from(del_umi)) as ArrayRef), // umi number (key)
        ("pos", Arc::new(UInt32Array::from(del_pos)) as ArrayRef), // 0-indexed position of first deleted reference base
        ("len", Arc::new(UInt32Array::from(del_len)) as ArrayRef), // deletion length
        ("obs", Arc::new(UInt64Array::from(del_obs)) as ArrayRef), // number of reads with this deletion
        ("hq", Arc::new(UInt64Array::from(del_hq)) as ArrayRef),   // number of reads that did NOT observe the deletion (high confidence)
        ("lq", Arc::new(UInt64Array::from(del_lq)) as ArrayRef),   // number of reads that did NOT observe the deletion (low confidence)
    ]).unwrap();
    let mut writer = ArrowWriter::try_new(del_file, del_batch.schema(), None).unwrap();
    writer.write(&del_batch).unwrap();
    writer.close().unwrap();
    
    eprintln!("Done!");
}
