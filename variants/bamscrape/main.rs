use std::fs::File;
use itertools::izip;
use std::str::from_utf8; //  return &str from &[u8]
// for bam
use rust_htslib::bam::{self,Read};
use rust_htslib::bam::record::Cigar; // enum
use rust_htslib::bam::record::CigarString; // struct Vec<Cigar>
// for .h5
use hdf5;

mod helpers;
use helpers::*;

type ReadNum = u32; // BAM line number
type SeqPos = usize;  // fasta byte position (must keep at usize)
type WLi = u32; // whitelist index (number of possible whitelist strings)

fn main() {
    // Read the BAM and FASTA path from the command line
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        panic!("Usage: bam_path fasta_path");
    }
    let bam_path = &args[1];
    let fasta_path = &args[2];
    
    // Load the BAM
    let mut bam = bam::Reader::from_path(&bam_path).expect("ERROR: BAM loading failed");
    // Get the list of RNAMEs from the BAM header
    let header = bam::Header::from_template(bam.header()).to_hashmap();
    let sq: Vec<String> = header.get("SQ").expect("ERROR: could not find @SQ in the header")
                                .iter().map(|e| e.get("SN").expect("SN key not found").clone())
                                .collect();

    // Load the FASTA as a memory-mapped file
    let fasta = File::open(&fasta_path).expect("ERROR: FASTA loading failed");
    let mmap = unsafe { memmap2::Mmap::map(&fasta).expect("ERROR: Mmap failed") };
    
    // Load the FASTA index file
    let index_path = [&fasta_path, ".fai"].join("");
    let map = load_fai(&index_path);
    
    // Summary of loaded data:
    // - bam is a rust_htslib::bam::Reader
    // - sq is a Vec<String> that maps tid -> RNAME
    // - map is a HashMap<String, usize> that maps RNAME -> fasta byte
    // - mmap is a Mmap containing the fasta
    assert!(sq.len() == map.len()); // techincally don't need
    assert!(mmap.len() <= SeqPos::MAX as usize);
    
    // Create string whitelists
    let mut cb_whitelist = Whitelist::new();
    let mut ub_whitelist = Whitelist::new();
    let mut ins_whitelist = Whitelist::new();
    let mut sc_whitelist = Whitelist::new();
    // Create data structures
    let mut data_read: Vec<ReadNum> = Vec::new();
    let mut data_flag: Vec<u16> = Vec::new();      // [2]
    let mut data_rname_i: Vec<usize> = Vec::new(); // [3]
    let mut data_pos: Vec<SeqPos> = Vec::new();    // [4]
    let mut data_mapq: Vec<u8> = Vec::new();       // [5]
    let mut data_cb: Vec<WLi> = Vec::new();        // CB
    let mut data_ub: Vec<WLi> = Vec::new();        // UB
    let mut data_xf: Vec<u8> = Vec::new();         // xf
    let mut data_re: Vec<u8> = Vec::new();         // RE
    let mut data_ts: Vec<SeqPos> = Vec::new();     // ts
    let mut data_pa: Vec<SeqPos> = Vec::new();     // pa
    // Create snv structures (M)
    let mut snv_read: Vec<ReadNum> = Vec::new();
    let mut snv_pos:  Vec<SeqPos> = Vec::new();
    let mut snv_ref: Vec<u8> = Vec::new();
    let mut snv_alt: Vec<u8> = Vec::new();
    let mut snv_qual: Vec<u8> = Vec::new();
    // Create insertion structures (I)
    let mut ins_read: Vec<ReadNum> = Vec::new();
    let mut ins_pos:  Vec<SeqPos> = Vec::new();
    let mut ins_str:  Vec<WLi> = Vec::new();
    // Create deletion structures (D)
    let mut del_read: Vec<ReadNum> = Vec::new();
    let mut del_pos:  Vec<SeqPos> = Vec::new();
    let mut del_len:  Vec<SeqPos> = Vec::new();
    // Create reference skip structures (N)
    let mut refskip_read: Vec<ReadNum> = Vec::new();
    let mut refskip_pos: Vec<SeqPos> = Vec::new();
    let mut refskip_len: Vec<SeqPos> = Vec::new();
    // Create match interval structures (M)
    let mut match_read: Vec<ReadNum> = Vec::new();
    let mut match_start: Vec<SeqPos> = Vec::new();
    let mut match_end: Vec<SeqPos> = Vec::new();
    // Create soft-clipped sequence structures (S)
    let mut sc_read: Vec<ReadNum> = Vec::new();
    let mut sc_pos:  Vec<SeqPos> = Vec::new();
    let mut sc_str:  Vec<WLi> = Vec::new();
    
    // Loop through the BAM records
    let mut read: ReadNum = 0;
    let mut no_cb: ReadNum = 0;
    let mut no_ub: ReadNum = 0;
    let mut no_rname: ReadNum = 0;
    let mut record = bam::Record::new();
    while let Some(r) = bam.read(&mut record) {
        r.expect("Failed to parse record");
        read = read.checked_add(1).expect("ERROR: BAM too long, must increase ReadNum");
        
        /* load the bam tags */
        
        // get the CB tag
        let cb: &str = match record.aux(b"CB") {
            Ok(value) => if let bam::record::Aux::String(s) = value { s } else {panic!("CB tag is not a String")},
            Err(_) => {no_cb+=1 ; ""}, // empty string, if CB tag does not exist
        };
        
        // get the UB (UMI) tag
        let ub: &str = match record.aux(b"UB") {
            Ok(value) => if let bam::record::Aux::String(s) = value { s } else {panic!("UB tag is not a String")},
            Err(_) => {no_ub+=1 ; ""}, // empty string, if UB tag does not exist
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
        
        // get the ts tag (number of trimmed TSO nucleotides)
        let ts: usize = match record.aux(b"ts") {
            Ok(value) => match value {
                bam::record::Aux::I32(s) => s as usize,
                bam::record::Aux::U8(s) => s as usize,
                _ => panic!("ts tag is not an I32 or U8"),
            },
            Err(_) => 0 as usize, // 0 trimmed by default
        };
        
        // get the pa tag (number of trimmed poly-A nucleotides)
        let pa: usize = match record.aux(b"pa") {
            Ok(value) => match value {
                bam::record::Aux::I32(s) => s as usize,
                bam::record::Aux::U8(s) => s as usize,
                _ => panic!("pa tag is not an I32 or U8"),
            },
            Err(_) => 0 as usize, // 0 trimmed by default
        };
        
        /* load the alignment section */
        
        // get the RNAME [3]
        let tid: i32 = record.tid();
        if tid < 0 { no_rname+=1 ; continue } // tid is -1 if no alignment
        if cb=="" || ub=="" { continue } // at this point, continue if there is not enough information to record the record
        let rname_i = tid as usize;
        let rname: &str = &sq[rname_i];
        let byte: usize = *map.get(rname).expect(&format!("ERROR: '{rname}' not found in fasta index file"));
        
        // get the FLAG [2], POS [4], and MAPQ [5]
        let flag: u16 = record.flags();
        let pos: SeqPos = record.pos().try_into().expect("ERROR: SeqPos is too small");
        let mapq: u8 = record.mapq();
        
        // write to data
        data_read.push(read); // ReadNum
        data_flag.push(flag); // u16 [2]
        data_rname_i.push(rname_i); // usize [3]
        data_pos.push(pos); // SeqPos [4]
        data_mapq.push(mapq); // u8 [5]
        data_cb.push(cb_whitelist.get(&cb)); // WLi CB
        data_ub.push(ub_whitelist.get(&ub)); // WLi UB
        data_xf.push(xf); // u8 xf
        data_re.push(re); // u8 RE
        data_ts.push(ts); // usize ts
        data_pa.push(pa); // usize pa
        
        // Process the CIGAR, SEQ, and QUAL for variants
        let cigarstring: CigarString = record.cigar().take(); // [6]
        let seq: Vec<u8> = record.seq().as_bytes(); // [10]
        let qual: &[u8] = record.qual(); // [11]
        let mut seq_i: usize = 0; // query
        let mut mmap_i: usize = byte + pos as usize; // reference
        for cigar in cigarstring.iter() {
            match cigar {
                Cigar::Match(v) => { // M
                    let n = *v as usize ;
                    let seq_iter = (&seq[seq_i..(seq_i + n)]).iter();
                    let ref_iter = (&mmap[mmap_i..(mmap_i + n)]).iter();
                    for (i, (&seq_byte, &ref_byte)) in izip!(seq_iter,ref_iter).enumerate() {
                        if seq_byte != ref_byte {
                            snv_read.push(read);
                            snv_pos.push(mmap_i + i - byte);
                            snv_ref.push(ref_byte);
                            snv_alt.push(seq_byte);
                            snv_qual.push(qual[seq_i + i]);
                        }
                    }
                    match_read.push(read);
                    match_start.push(mmap_i - byte);
                    match_end.push(mmap_i + n - byte);
                    // println!("{:?}", cigarstring) ;
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
                Cigar::SoftClip(v) => { // S
                    let n = *v as usize;
                    sc_read.push(read);
                    sc_pos.push(mmap_i - byte);
                    sc_str.push(sc_whitelist.get(from_utf8(&seq[seq_i..(seq_i + n)]).expect("E")));
                    seq_i += n;
                },
                Cigar::HardClip(_) => panic!("Not implemented cigar operator 'H'"), // H
                Cigar::Pad(_) => panic!("Not implemented cigar operator 'P'"),      // P
                Cigar::Equal(_) => panic!("Not implemented cigar operator '='"),    // =
                Cigar::Diff(_) => panic!("Not implemented cigar operator 'X'"),     // X
            }
        }
    }
    
    // panic!("skipping writing output");
    
    // write to an .h5 file
    let file = hdf5::File::create("reads.h5").expect("Could not create .h5 output file");
    // write string whitelists
    let whitelist_group = file.create_group("whitelist").expect(".h5 error");
    cb_whitelist.into_vector().write_vector(&whitelist_group, "cb");   // list of all observed cell barcodes
    ub_whitelist.into_vector().write_vector(&whitelist_group, "ub");   // list of all observed UMI barcodes
    ins_whitelist.into_vector().write_vector(&whitelist_group, "ins"); // list of all inserted strings
    sc_whitelist.into_vector().write_vector(&whitelist_group, "sc");   // list of all soft-clipped sequences
    sq.write_vector(&whitelist_group,"rname");                         // list of all RNAME
    // write general data
    let data_group = file.create_group("data").expect(".h5 error");
    assert_all_same(&[data_read.len(), data_flag.len(), data_rname_i.len(), data_pos.len(), data_mapq.len(), data_cb.len(), data_ub.len(),data_xf.len(),data_re.len(),data_ts.len(),data_pa.len()]);
    data_read.write_vector(&data_group, "read");       // bam line number (1-indexed)
    data_flag.write_vector(&data_group, "flag");       // FLAG
    data_rname_i.write_vector(&data_group, "rname_i"); // RNAME index
    data_pos.write_vector(&data_group, "pos");         // POS
    data_mapq.write_vector(&data_group, "mapq");       // MAPQ
    data_cb.write_vector(&data_group, "cb_i");         // cb whitelist index (0-indexed)
    data_ub.write_vector(&data_group, "ub_i");         // ub whitelist index (0-indexed)
    data_xf.write_vector(&data_group, "xf");           // xf tag
    data_re.write_vector(&data_group, "re");           // RE tag
    data_ts.write_vector(&data_group, "ts");           // ts tag
    data_pa.write_vector(&data_group, "pa");           // pa tag
    // write snv data
    let snv_group = file.create_group("snv").expect(".h5 error");
    assert_all_same(&[snv_read.len(), snv_pos.len(), snv_ref.len(), snv_alt.len(), snv_qual.len()]);
    snv_read.write_vector(&snv_group, "read");  // bam line number (1-indexed)
    snv_pos.write_vector(&snv_group, "pos");    // 0-indexed reference position of snv
    snv_ref.write_vector(&snv_group, "ref");    // reference base
    snv_alt.write_vector(&snv_group, "alt");    // alternate base
    snv_qual.write_vector(&snv_group, "qual");  // base quality
    // write insertion data
    let ins_group = file.create_group("ins").expect(".h5 error");
    assert_all_same(&[ins_read.len(), ins_pos.len(), ins_str.len()]);
    ins_read.write_vector(&ins_group, "read");            // bam line number (1-indexed)
    ins_pos.write_vector(&ins_group, "pos");              // 0-indexed position of reference base after insertion
    ins_str.write_vector(&ins_group, "str_i");            // index into whitelist/ins (0-indexed)
    // write deletion data
    let del_group = file.create_group("del").expect(".h5 error");
    assert_all_same(&[del_read.len(), del_pos.len(), del_len.len()]);
    del_read.write_vector(&del_group, "read");            // bam line number (1-indexed)
    del_pos.write_vector(&del_group, "pos");              // 0-indexed position of first deleted reference base
    del_len.write_vector(&del_group, "len");              // number of deleted bases
    // write refskip (N) data
    let refskip_group = file.create_group("refskip").expect(".h5 error");
    assert_all_same(&[refskip_read.len(), refskip_pos.len(), refskip_len.len()]);
    refskip_read.write_vector(&refskip_group, "read");    // bam line number (1-indexed)
    refskip_pos.write_vector(&refskip_group, "pos");      // 0-indexed position of first skipped reference base
    refskip_len.write_vector(&refskip_group, "len");      // number of skipped bases
    // write matching interval (M) data
    let match_group = file.create_group("match").expect(".h5 error");
    assert_all_same(&[match_read.len(), match_start.len(), match_end.len()]);
    match_read.write_vector(&match_group, "read");        // bam line number (1-indexed)
    match_start.write_vector(&match_group, "start");      // 0-indexed position of range start (inclusive)
    match_end.write_vector(&match_group, "end");          // 0-indexed position of range end (exclusive)
    // write soft-clipped (S) data
    let sc_group = file.create_group("softclip").expect(".h5 error");
    assert_all_same(&[sc_read.len(), sc_pos.len(), sc_str.len()]);
    sc_read.write_vector(&sc_group, "read");              // bam line number (1-indexed)
    sc_pos.write_vector(&sc_group, "start");              // 0-indexed position of reference base after skip
    sc_str.write_vector(&sc_group, "end");                // index into whitelist/sc (0-indexed)
    // write metadata
    let meta_group = file.create_group("metadata").expect(".h5 error");
    let names= vec!["reads","no_cb","no_ub","no_rname"];
    let counts = vec![read, no_cb, no_ub, no_rname];
    assert_all_same(&[names.len(), counts.len()]);
    names.write_vector(&meta_group,"names");
    counts.write_vector(&meta_group, "read");
    
    println!("DONE");
}
