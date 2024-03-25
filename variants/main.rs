use std::env; // for CLA
use std::path::Path;
use std::fs::File;
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::str::from_utf8; //  return &str
// for bam
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar; // enum
use rust_htslib::bam::record::CigarString; // struct Vec<Cigar>
// for Mmap
use memmap2::Mmap;
// for .h5
use hdf5;
use hdf5::types::VarLenAscii;

struct Whitelist {
    map: HashMap<String, usize>
}
impl Whitelist {
    // create an empty whitelist
    fn new() -> Whitelist { Whitelist { map: HashMap::new() } }
    // return the index of the string in the whitelist (add if needed)
    fn get(&mut self, string: &str) -> usize {
        match self.map.get(string) {
            Some(&val) => val,
            None => { let n = self.map.len(); self.map.insert(string.to_string(), n); n }
        }
    }
    // return a vector of strings (in order)
    fn into_vector(self) -> Vec<String> {
        let mut pairs = self.map.into_iter().collect::<Vec<(String, usize)>>();
        pairs.sort_by(|a, b| a.1.cmp(&b.1));
        pairs.into_iter().map(|(key, _)| key).collect::<Vec<String>>()
    }
    // return a vector of strings (in order), able to be added to an .h5 file
    fn into_H5Type(self) -> Vec<VarLenAscii> {
      self.into_vector().into_iter().map(|s| VarLenAscii::from_ascii(&s).expect("ERROR: non-ascii character")).collect()
    }
}

// create a (RNAME -> fasta byte offset) map from the .fai file
fn load_fai(index_path: &str) -> HashMap<String, usize> {
    assert!(Path::new(&index_path).exists(), "ERROR: fasta index file not found");
    let mut map: HashMap<String, usize> = HashMap::new();
    let mut index_file = ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(File::open(&index_path).expect("ERROR: could not open index file"));
    for result in index_file.records() {
        let record = result.expect("ERROR: could not parse index file");
        if let (Some(key), Some(value)) = (record.get(0), record.get(2)) { // col 1 is RNAME, col 3 is fasta byte offset
            map.insert(key.to_string(), value.parse::<usize>().expect("ERROR: .fai column 3 could not be read as usize"));
        } else { panic!("ERROR: .fai file missing column data"); }
    }
    return map;
}

// write a vector to an .h5 group
enum DownsizedVector {
    U8Vec(Vec<u8>),
    U16Vec(Vec<u16>),
    U32Vec(Vec<u32>),
    U64Vec(Vec<u64>),
}
fn downsize(vec: Vec<usize>) -> DownsizedVector {
    let max: usize = *vec.iter().max().unwrap_or(&0);
    match max {
        0..=255            => DownsizedVector::U8Vec(vec.into_iter().map(|x| x as u8).collect()),
        256..=65535        => DownsizedVector::U16Vec(vec.into_iter().map(|x| x as u16).collect()),
        65536..=4294967295 => DownsizedVector::U32Vec(vec.into_iter().map(|x| x as u32).collect()),
        _                  => DownsizedVector::U64Vec(vec.into_iter().map(|x| x as u64).collect()),
    }
}
fn write_vector(group: &hdf5::Group, vec: DownsizedVector, name: &str) {
    match vec {
        DownsizedVector::U8Vec(v) => group.new_dataset_builder().with_data(&v).create(name).expect(".h5 error"),
        DownsizedVector::U16Vec(v) => group.new_dataset_builder().with_data(&v).create(name).expect(".h5 error"),
        DownsizedVector::U32Vec(v) => group.new_dataset_builder().with_data(&v).create(name).expect(".h5 error"),
        DownsizedVector::U64Vec(v) => group.new_dataset_builder().with_data(&v).create(name).expect(".h5 error"),
    };
}

fn main() {
    // Read the BAM and FASTA path from the command line
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        println!("Usage: bam_path fasta_path");
        // assert!(false);
    }
    // let bam_path = &args[1];
    // let fasta_path = &args[2];
    
    let bam_path = "/home/nsachdev/analyses/0data/possorted_genome_bam.bam";
    // let bam_path = "/home/nsachdev/analyses/8snpscrape/little.bam";
    let fasta_path = "/home/nsachdev/analyses/8snpscrape/removefastawhitespace/fasta_nowhitespace/refdata-gex-GRCh38-2020-A.fa";

    // Load the FASTA as a memory-mapped file
    assert!(Path::new(&fasta_path).exists(), "ERROR: input fasta file not found");
    let fasta = File::open(&fasta_path).expect("ERROR: FASTA loading failed");
    let mmap = unsafe { Mmap::map(&fasta).expect("ERROR: Mmap failed") };
    
    // Load the FASTA index file
    let index_path = [&fasta_path, ".fai"].join("");
    let map = load_fai(&index_path);
    
    // Load the BAM
    let mut bam = bam::Reader::from_path(&bam_path).expect("ERROR: BAM loading failed");
    
    // Get the list of RNAMEs from the BAM header
    let header = bam::Header::from_template(bam.header()).to_hashmap();
    let sq = header.get("SQ").expect("ERROR: could not find @SQ in the header")
                 .iter().map(|e| e.get("SN").expect("SN key not found").clone())
                 .collect::<Vec<String>>();
    
    // Summary of loaded data:
    // - bam is a rust_htslib::bam::Reader
    // - sq is a Vec<String> that maps tid -> RNAME
    // - map is a HashMap<String, usize> that maps RNAME -> fasta byte
    // - mmap is a Mmap containing the fasta
    
    // Create string whitelists
    let mut cb_whitelist = Whitelist::new();
    let mut ub_whitelist = Whitelist::new();
    let mut ins_whitelist = Whitelist::new();
    // Create data structures
    let mut data_read: Vec<usize> = Vec::new();
    let mut data_cb: Vec<usize> = Vec::new();
    let mut data_ub: Vec<usize> = Vec::new();
    let mut data_xf: Vec<u8> = Vec::new();
    let mut data_re: Vec<u8> = Vec::new();
    let mut data_flag: Vec<u16> = Vec::new();  // [2]
    let mut data_rname: Vec<u32> = Vec::new(); // [3]
    let mut data_pos: Vec<usize> = Vec::new(); // [4]
    let mut data_mapq: Vec<u8> = Vec::new();   // [5]
    // Create snv structures
    let mut snv_read: Vec<usize> = Vec::new();
    let mut snv_pos:  Vec<usize> = Vec::new();
    let mut snv_ref: Vec<u8> = Vec::new();
    let mut snv_alt: Vec<u8> = Vec::new();
    let mut snv_qual: Vec<u8> = Vec::new();
    // Create insertion structures
    let mut ins_read: Vec<usize> = Vec::new();
    let mut ins_pos:  Vec<usize> = Vec::new();
    let mut ins_str:  Vec<usize> = Vec::new();
    // Create deletion structures
    let mut del_read: Vec<usize> = Vec::new();
    let mut del_pos:  Vec<usize> = Vec::new();
    let mut del_len:  Vec<usize> = Vec::new();
    // Create reference skip structures
    let mut refskip_read: Vec<usize> = Vec::new();
    let mut refskip_pos: Vec<usize> = Vec::new();
    let mut refskip_len: Vec<usize> = Vec::new();
    // Create match interval structures
    let mut match_read: Vec<usize> = Vec::new();
    let mut match_start: Vec<usize> = Vec::new();
    let mut match_end: Vec<usize> = Vec::new();
    
    // Loop through the BAM records
    let mut read: usize = 0;
    let mut record = bam::Record::new();
    while let Some(r) = bam.read(&mut record) {
        r.expect("Failed to parse record");
        read += 1;
        
        // get the CB tag
        let cb: &str = match record.aux(b"CB") {
            Ok(value) => if let bam::record::Aux::String(s) = value { s } else {panic!("CB tag is not a String")},
            Err(_) => continue, // if the CB tag does not exist, continue to next loop iteration
        };
        
        // get the UB (UMI) tag
        let ub: &str = match record.aux(b"UB") {
            Ok(value) => if let bam::record::Aux::String(s) = value { s } else {panic!("UB tag is not a String")},
            Err(_) => continue, // if the UB tag does not exist, continue to next loop iteration
        };
        
        // get the xf tag (extra alignment flags)
        let xf: u8 = match record.aux(b"xf") {
            Ok(value) => match value {
                bam::record::Aux::I32(s) => s as u8,
                bam::record::Aux::U8(s) => s,
                _ => panic!("xf tag is not an I32 or U8"),
            },
            Err(_) => 255u8, // xf only goes up to 63, so assign 255 if blank
        };
        
        // get the RE tag (alignment region type)
        let re: u8 = match record.aux(b"RE") {
            Ok(value) => if let bam::record::Aux::Char(s) = value { s } else {panic!("RE tag is not a Char")},
            Err(_) => 'X' as u8, // RE is either E/N/I, so assign X if blank
        };
        
        // get the RNAME [3]
        let tid: i32 = record.tid();
        if tid < 0 { continue } // tid is -1 if no alignment
        let rname: &str = &sq[tid as usize];
        let byte: usize = *map.get(rname).expect(&format!("ERROR: '{rname}' not found in fasta index file"));
        
        // get the FLAG [2], POS [4], and MAPQ [5]
        let flag: u16 = record.flags();
        let pos: usize = record.pos() as usize;
        let mapq: u8 = record.mapq();
        
        // write to data
        data_read.push(read);
        data_cb.push(cb_whitelist.get(&cb));
        data_ub.push(ub_whitelist.get(&ub));
        data_xf.push(xf);
        data_re.push(re);
        data_flag.push(flag);
        data_rname.push(tid as u32);
        data_pos.push(pos);
        data_mapq.push(mapq);
        
        // Process the CIGAR, SEQ, and QUAL for variants
        let cigarstring: CigarString = record.cigar().take(); // [6]
        let seq: Vec<u8> = record.seq().as_bytes(); // [10]
        let qual: &[u8] = record.qual(); // [11]
        let mut seq_i: usize = 0; // query
        let mut mmap_i: usize = byte + pos; // reference
        for cigar in cigarstring.iter() {
            match cigar {
                Cigar::Match(v) => { // M
                    let n = *v as usize ;
                    let seq_iter = (&seq[seq_i..(seq_i + n)]).iter();
                    let ref_iter = (&mmap[mmap_i..(mmap_i + n)]).iter();
                    for (i, (&byte_seq, &byte_ref)) in seq_iter.zip(ref_iter).enumerate() {
                        if byte_seq != byte_ref {
                            snv_read.push(read);
                            snv_pos.push(mmap_i + i - byte);
                            snv_ref.push(byte_ref);
                            snv_alt.push(byte_seq);
                            snv_qual.push(qual[seq_i + i]);
                        }
                    }
                    match_read.push(read);
                    match_start.push(mmap_i - byte);
                    match_end.push(mmap_i + n - byte);
                    // println!("{:?}", from_utf8(&seq[seq_i..(seq_i + n)]).expect("E") ) ;
                    // println!("{:?}", from_utf8(&mmap[mmap_i..(mmap_i + n)]).expect("E") ) ;
                    seq_i += n;
                    mmap_i += n;
                },
                Cigar::Ins(v) => { // I
                    let n = *v as usize;
                    ins_read.push(read);
                    ins_pos.push(mmap_i - byte);
                    ins_str.push(ins_whitelist.get(from_utf8(&seq[seq_i..(seq_i + n)]).expect("E")));
                    seq_i += n;
                }, 
                Cigar::Del(v) => { // D
                    let n = *v as usize;
                    del_read.push(read);
                    del_pos.push(mmap_i - byte);
                    del_len.push(n);
                    mmap_i += n;
                },
                Cigar::RefSkip(v) => { // N
                    let n = *v as usize;
                    refskip_read.push(read);
                    refskip_pos.push(mmap_i - byte);
                    refskip_len.push(n);
                    mmap_i += n;
                }, 
                Cigar::SoftClip(v) => seq_i += *v as usize, // S
                Cigar::HardClip(_) => panic!("Not implemented cigar operator 'H'"), // H
                Cigar::Pad(_) => panic!("Not implemented cigar operator 'P'"),      // P
                Cigar::Equal(_) => panic!("Not implemented cigar operator '='"),    // =
                Cigar::Diff(_) => panic!("Not implemented cigar operator 'X'"),     // X
            }
        }
    }
    
    //assert!(false);
    
    // write to an .h5 file
    let file = hdf5::File::create("raw.h5").expect("Could not create .h5 output file");
    
    // write string whitelists
    let whitelist_group = file.create_group("whitelist").expect(".h5 error");
    whitelist_group.new_dataset_builder().with_data(&cb_whitelist.into_H5Type()).create("cb").expect(".h5 error"); // list of all observed cell barcodes
    whitelist_group.new_dataset_builder().with_data(&ub_whitelist.into_H5Type()).create("ub").expect(".h5 error"); // list of all observed UMI barcodes
    whitelist_group.new_dataset_builder().with_data(&ins_whitelist.into_H5Type()).create("ins").expect(".h5 error"); // list of all inserted strings
    let sq_h5: Vec<VarLenAscii> = sq.into_iter().map(|s| VarLenAscii::from_ascii(&s).expect("ERROR: non-ascii character")).collect();
    whitelist_group.new_dataset_builder().with_data(&sq_h5).create("rname").expect(".h5 error"); // list of all RNAME
    
    // write general data
    let data_group = file.create_group("data").expect(".h5 error");
    assert!(data_read.len() == data_cb.len()); assert!(data_cb.len() == data_ub.len()); assert!(data_ub.len() == data_xf.len());
    assert!(data_xf.len() == data_re.len()); assert!(data_re.len() == data_flag.len());
    assert!(data_flag.len() == data_rname.len()); assert!(data_rname.len() == data_mapq.len());
    write_vector(&data_group, downsize(data_read), "read");                    // bam line number (1-indexed)
    write_vector(&data_group, downsize(data_cb), "cb_i");                      // cb whitelist index (0-indexed)
    write_vector(&data_group, downsize(data_ub), "ub_i");                      // ub whitelist index (0-indexed)
    write_vector(&data_group, DownsizedVector::U8Vec(data_xf), "xf");          // xf tag
    write_vector(&data_group, DownsizedVector::U8Vec(data_re), "re");          // RE tag
    write_vector(&data_group, DownsizedVector::U16Vec(data_flag), "flag");     // FLAG
    write_vector(&data_group, DownsizedVector::U32Vec(data_rname), "rname_i"); // RNAME index
    write_vector(&data_group, downsize(data_pos), "pos");                      // POS
    write_vector(&data_group, DownsizedVector::U8Vec(data_mapq), "mapq");      // MAPQ
    
    // write snv data
    let snv_group = file.create_group("snv").expect(".h5 error");
    assert!(snv_read.len() == snv_pos.len()); assert!(snv_pos.len() == snv_ref.len());
    assert!(snv_ref.len() == snv_alt.len()); assert!(snv_alt.len() == snv_qual.len());
    write_vector(&snv_group, downsize(snv_read), "read");               // bam line number (1-indexed)
    write_vector(&snv_group, downsize(snv_pos), "pos");                 // 0-indexed reference position of snv
    write_vector(&snv_group, DownsizedVector::U8Vec(snv_ref), "ref");   // reference base
    write_vector(&snv_group, DownsizedVector::U8Vec(snv_alt), "alt");   // alternate base
    write_vector(&snv_group, DownsizedVector::U8Vec(snv_qual), "qual"); // base quality
    
    // write insertion data
    let ins_group = file.create_group("ins").expect(".h5 error");
    assert!(ins_read.len() == ins_pos.len()); assert!(ins_pos.len() == ins_str.len());
    write_vector(&ins_group, downsize(ins_read), "read"); // bam line number (1-indexed)
    write_vector(&ins_group, downsize(ins_pos), "pos");   // 0-indexed position of reference base after insertion
    write_vector(&ins_group, downsize(ins_str), "str_i"); // index into whitelist/ins (0-indexed)
    
    // write deletion data
    let del_group = file.create_group("del").expect(".h5 error");
    assert!(del_read.len() == del_pos.len()); assert!(del_pos.len() == del_len.len());
    write_vector(&del_group, downsize(del_read), "read"); // bam line number (1-indexed)
    write_vector(&del_group, downsize(del_pos), "pos");   // 0-indexed position of first deleted reference base
    write_vector(&del_group, downsize(del_len), "len");   // number of deleted bases
    
    // write refskip (N) data
    let refskip_group = file.create_group("refskip").expect(".h5 error");
    assert!(refskip_read.len() == refskip_pos.len()); assert!(refskip_pos.len() == refskip_len.len());
    write_vector(&refskip_group, downsize(refskip_read), "read"); // bam line number (1-indexed)
    write_vector(&refskip_group, downsize(refskip_pos), "pos");   // 0-indexed position of first skipped reference base
    write_vector(&refskip_group, downsize(refskip_len), "len");   // number of skipped bases
    
    // write matching interval (M) data
    let match_group = file.create_group("match").expect(".h5 error");
    assert!(match_read.len() == match_start.len()); assert!(match_start.len() == match_end.len());
    write_vector(&match_group, downsize(match_read), "read");   // bam line number (1-indexed)
    write_vector(&match_group, downsize(match_start), "start"); // 0-indexed position of range start (inclusive)
    write_vector(&match_group, downsize(match_end), "end");     // 0-indexed position of range end (exclusive)
    
    println!("DONE");
}
