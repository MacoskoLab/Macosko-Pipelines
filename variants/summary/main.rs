use std::collections::HashMap;
// use hdf5::types::VarLenAscii;
use itertools::izip;
use hdf5;
// for writing output
use csv::{WriterBuilder, Terminator};
use flate2::{write::GzEncoder, Compression};
use std::fs::File;

const QUAL_CUTOFF: u8 = 30; // 2 12 26 34

#[derive(Debug)]
struct UMIinfo {
    // snv
    pos: Vec<u32>,
    alt: Vec<u8>,
    qual: Vec<u8>,
    // match
    start: Vec<u32>,
    end: Vec<u32>,
}
impl UMIinfo {
    fn new() -> UMIinfo { UMIinfo { pos: Vec::new(), alt: Vec::new(), qual: Vec::new(), start: Vec::new(), end:Vec::new() } }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        panic!("Usage: reads.h5_path");
    }
    let file = hdf5::File::open(&args[1]).expect("Could not open reads.h5 input file");
    
    // let cb_whitelist: Vec<VarLenAscii> = file.dataset("whitelist/cb").expect("not found").read_raw().expect("E");
    // let ub_whitelist: Vec<VarLenAscii> = file.dataset("whitelist/ub").expect("not found").read_raw().expect("E");
    // let ins_whitelist: Vec<VarLenAscii> = file.dataset("whitelist/ins").expect("not found").read_raw().expect("E");
    // let sc_whitelist: Vec<VarLenAscii> = file.dataset("whitelist/sc").expect("not found").read_raw().expect("E");
    // let rname_whitelist: Vec<VarLenAscii> = file.dataset("whitelist/rname").expect("not found").read_raw().expect("E");
    
    let data_read: Vec<u32> = file.dataset("data/read").expect("not found").read_raw().expect("E");
    // let data_flag: Vec<u32> = file.dataset("data/flag").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_flag.len());
    let data_rname_i: Vec<u32> = file.dataset("data/rname_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_rname_i.len());
    // let data_pos: Vec<u64> = file.dataset("data/pos").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_pos.len());
    // let data_mapq: Vec<u32> = file.dataset("data/mapq").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_mapq.len());
    let data_cb_i: Vec<u32> = file.dataset("data/cb_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_cb_i.len());
    let data_ub_i: Vec<u32> = file.dataset("data/ub_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_ub_i.len());
    // let data_re: Vec<u64> = file.dataset("data/re").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_re.len());
    // let data_xf: Vec<u64> = file.dataset("data/xf").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_xf.len());
    // let data_ts: Vec<u64> = file.dataset("data/ts").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_ts.len());
    // let data_pa: Vec<u64> = file.dataset("data/pa").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_pa.len());
    
    let snv_read: Vec<u32> = file.dataset("snv/read").expect("not found").read_raw().expect("E");
    let snv_pos: Vec<u32> = file.dataset("snv/pos").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_pos.len());
    // let snv_ref: Vec<u8> = file.dataset("snv/ref").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_ref.len());
    let snv_alt: Vec<u8> = file.dataset("snv/alt").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_alt.len());
    let snv_qual: Vec<u8> = file.dataset("snv/qual").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_qual.len());
    
    // let ins_read: Vec<u64> = file.dataset("ins/read").expect("not found").read_raw().expect("E");
    // let ins_pos: Vec<u64> = file.dataset("ins/pos").expect("not found").read_raw().expect("E"); assert!(ins_read.len() == ins_pos.len());
    // let ins_str_i: Vec<u64> = file.dataset("ins/str_i").expect("not found").read_raw().expect("E"); assert!(ins_read.len() == ins_str_i.len());
    
    // let del_read: Vec<u64> = file.dataset("del/read").expect("not found").read_raw().expect("E");
    // let del_pos: Vec<u64> = file.dataset("del/pos").expect("not found").read_raw().expect("E"); assert!(del_read.len() == del_pos.len());
    // let del_len: Vec<u64> = file.dataset("del/len").expect("not found").read_raw().expect("E"); assert!(del_read.len() == del_len.len());
    
    // let refskip_read: Vec<u64> = file.dataset("refskip/read").expect("not found").read_raw().expect("E");
    // let refskip_pos: Vec<u64> = file.dataset("refskip/pos").expect("not found").read_raw().expect("E"); assert!(refskip_read.len() == refskip_pos.len());
    // let refskip_len: Vec<u64> = file.dataset("refskip/len").expect("not found").read_raw().expect("E"); assert!(refskip_read.len() == refskip_len.len());
    
    let match_read: Vec<u32> = file.dataset("match/read").expect("not found").read_raw().expect("E");
    let match_start: Vec<u32> = file.dataset("match/start").expect("not found").read_raw().expect("E"); assert!(match_read.len() == match_start.len());
    let match_end: Vec<u32> = file.dataset("match/end").expect("not found").read_raw().expect("E"); assert!(match_read.len() == match_end.len());
    
    // for trans-splicing
    // let sc_read: Vec<u64> = file.dataset("softclip/read").expect("not found").read_raw().expect("E");
    // let sc_pos: Vec<u64> = file.dataset("softclip/pos").expect("not found").read_raw().expect("E"); assert!(sc_read.len() == sc_pos.len());
    // let sc_str_i: Vec<u64> = file.dataset("softclip/str_i").expect("not found").read_raw().expect("E"); assert!(sc_read.len() == sc_str_i.len());
    // let data_ts: Vec<u64> = file.dataset("data/ts").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_ts.len());
    // let data_pa: Vec<u64> = file.dataset("data/pa").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_pa.len());
    
    // load read2data_i and umi_table
    let mut read2data_i: HashMap<u32, usize> = HashMap::new(); // read -> data_i
    let mut umi_table: HashMap<(u32,u32,u32), UMIinfo> = HashMap::new(); // (cb_i, ub_i, rname_i) -> UMIinfo
    for i in 0..data_read.len() {
      read2data_i.insert(data_read[i], i);
      umi_table.entry((data_cb_i[i], data_ub_i[i], data_rname_i[i])).or_insert_with(|| UMIinfo::new());
    }
    
    // get snv per umi
    for i in 0..snv_read.len() {
        let data_i: usize = *read2data_i.get(&snv_read[i]).expect("E"); // get the index in "data"
        let key = (data_cb_i[data_i], data_ub_i[data_i], data_rname_i[data_i]); // get the umi
        let umiinfo: &mut UMIinfo = umi_table.get_mut(&key).expect("E");
        umiinfo.pos.push(snv_pos[i]);
        umiinfo.alt.push(snv_alt[i]);
        umiinfo.qual.push(snv_qual[i]);
    }
    
    // get match ranges per umi
    for i in 0..match_read.len() {
        let data_i: usize = *read2data_i.get(&match_read[i]).expect("E"); // get the index in "data"
        let key = (data_cb_i[data_i], data_ub_i[data_i], data_rname_i[data_i]); // get the umi
        let umiinfo: &mut UMIinfo = umi_table.get_mut(&key).expect("E");
        umiinfo.start.push(match_start[i]);
        umiinfo.end.push(match_end[i]);
    }
    
    // aggregate snv/umi counts for/against
    //              cb_i ub_i rname_i pos alt  hq   lq  total
    let mut ret: Vec<(u32, u32, u32, u32, u8, u32, u32, u32)> = Vec::new();
    for ((cb_i,umi_i,rname_i), umiinfo) in umi_table.iter() {
        // get the number of hq and lq reads for each snv
        let mut pos: Vec<u32> = Vec::new();
        let mut alt: Vec<u8> = Vec::new();
        let mut hq: Vec<u32> = Vec::new();
        let mut lq: Vec<u32> = Vec::new();
        let mut counts_map: HashMap<(u32, u8), (u32, u32)> = HashMap::new(); // (pos, alt) -> (hq, lq)
        for (&p, &a, &q) in izip!(&umiinfo.pos, &umiinfo.alt, &umiinfo.qual) {
            counts_map.entry((p, a))
                .and_modify(|e| if q > QUAL_CUTOFF { e.0 += 1 } else { e.1 += 1 })
                .or_insert(if q > QUAL_CUTOFF { (1, 0) } else { (0, 1) });
        }
        for ((p, a), (h, l)) in counts_map {
            pos.push(p); alt.push(a); hq.push(h); lq.push(l);
        }
        
        // get the total number of observations for each snv
        let total: Vec<u32> = pos.iter().map(|&p| {
            izip!(&umiinfo.start, &umiinfo.end).filter(|(&s, &e)| p >= s && p < e).count() as u32
        }).collect();
        
        // write to memory
        assert!(pos.len() == alt.len()); assert!(pos.len() == hq.len());
        assert!(pos.len() == lq.len()); assert!(pos.len() == total.len());
        for i in 0..pos.len() {
            ret.push((*cb_i, *umi_i, *rname_i, pos[i], alt[i], hq[i], lq[i], total[i]));
        }
    }
    
    // write the output
    println!("Writing output");
    let file = File::create("snv.tsv.gz").expect("ERROR: could not create snv.tsv.gz output file");
    let encoder = GzEncoder::new(file, Compression::default());
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::Any(b'\n'))
        .from_writer(encoder);
    writer.write_record(&["cb_i", "ub_i", "rname_i", "pos", "alt", "hq", "lq", "total"]).expect("E");
    for row in ret {
        writer.serialize(row).expect("E");
    }
    writer.flush().expect("E");
    
    println!("DONE");
}
