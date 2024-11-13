use std::collections::HashMap;
use std::collections::HashSet;
use num_traits::int::PrimInt;
use itertools::izip;

pub type ReadNum = u64;
pub type SeqPos = u32;
pub type WLi = u32;
pub type RNAMEi = u16;
pub type SNVi = usize;

// Data structure that aggregates all the information of a UMI
pub struct UMIinfo {
    pub reads: ReadNum,
    pub mapq: Vec<u8>,
    pub snv_i: HashMap<(SNVi, SeqPos), (ReadNum, ReadNum)>, // (snv_i, pos) -> (hq, lq)
    pub matches: Vec<(SeqPos, SeqPos)>, // (start, end)
}
impl UMIinfo {
    pub fn new() -> Self {
        Self { reads: 0, mapq: Vec::new(), snv_i: HashMap::new(), matches: Vec::new() }
    }
    
    // collapse match intervals into a consensus range
    pub fn merge_intervals(&self) -> Vec<(SeqPos, SeqPos)> {
        let mut intervals = self.matches.clone();
        intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged_intervals: Vec<(SeqPos, SeqPos)> = Vec::new();
        for interval in intervals {
            if let Some(last) = merged_intervals.last_mut() {
                if interval.0 <= last.1 {
                    last.1 = last.1.max(interval.1); // extend the current interval
                } else {
                    merged_intervals.push(interval); // no overlap, add a new interval
                }
            } else {
                merged_intervals.push(interval); // first interval
            }
        }
        return merged_intervals
    }
    
    // get the average value of a list of u8 (used for computing average mapq)
    pub fn mapq_avg(&self) -> u8 {
        let sum: u64 = self.mapq.iter().map(|&val| u64::from(val)).sum();
        let len: u64 = self.mapq.len() as u64;
        let avg = (sum + len/2) / len;
        return avg as u8
    }
}

// Data structure that enables storing "indexes into a list of SNV" instead of "SNV"
pub struct SNVTable {
    map: HashMap<(RNAMEi, SeqPos, u8), SNVi> // (rname_i, pos, alt) -> snv_i
}
impl SNVTable {
    // create an empty whitelist
    pub fn new() -> Self {
        Self { map: HashMap::new() }
    }
    // return the index of the string in the whitelist (add if needed)
    pub fn get(&mut self, snv: &(RNAMEi, SeqPos, u8)) -> SNVi {
        match self.map.get(snv) {
            Some(&val) => val,
            None => {
                let n: SNVi = self.map.len();
                self.map.insert(*snv, n);
                n
            }
        }
    }
    // return a tuple of vectors (in order)
    pub fn into_vectors(self) -> (Vec<SNVi>, Vec<RNAMEi>, Vec<SeqPos>, Vec<u8>) {
        let mut pairs: Vec<((RNAMEi, SeqPos, u8), SNVi)> = self.map.into_iter().collect();
        pairs.sort_by(|a, b| a.1.cmp(&b.1));
        let mut vec_snv_i: Vec<SNVi> = Vec::new();
        let mut vec_rname_i: Vec<RNAMEi> = Vec::new();
        let mut vec_pos: Vec<SeqPos> = Vec::new();
        let mut vec_alt: Vec<u8> = Vec::new();
        for ((rname_i, pos, alt), snv_i) in pairs {
            vec_snv_i.push(snv_i);
            vec_rname_i.push(rname_i);
            vec_pos.push(pos);
            vec_alt.push(alt);
        }
        return (vec_snv_i, vec_rname_i, vec_pos, vec_alt)
    }
}

// extract the strand from the flag
pub fn flag2strand(flag: u16) -> u8 {
  ((flag >> 4) & 1) as u8
}

// helper method to assert that all elements of a vector are the same
pub fn assert_all_same(vec: &[usize]) {
    if let Some(first) = vec.first() {
        for element in vec.iter() {
            assert_eq!(element, first, "ERROR: not all lengths are the same");
        }
    }
}

/* DATA */
pub struct DataEntry {
    pub read: ReadNum,
    pub flag: u16,
    pub rname_i: RNAMEi,
    pub mapq: u8,
    pub cb_i: WLi,
    pub ub_i: WLi,
}
pub fn load_data(file: &hdf5::File) -> Vec<DataEntry> {
    assert!(file.dataset("data/read").expect("not found").dtype().expect(".h5 error").size() <= size_of::<ReadNum>(), "ERROR: ReadNum too small");
    assert!(file.dataset("data/flag").expect("not found").dtype().expect(".h5 error").size() <= size_of::<u16>(), "ERROR: u16 too small");
    assert!(file.dataset("data/rname_i").expect("not found").dtype().expect(".h5 error").size() <= size_of::<RNAMEi>(), "ERROR: RNAMEi too small");
    assert!(file.dataset("data/mapq").expect("not found").dtype().expect(".h5 error").size() <= size_of::<u8>(), "ERROR: u8 too small");
    assert!(file.dataset("data/cb_i").expect("not found").dtype().expect(".h5 error").size() <= size_of::<WLi>(), "ERROR: WLi too small");
    assert!(file.dataset("data/ub_i").expect("not found").dtype().expect(".h5 error").size() <= size_of::<WLi>(), "ERROR: WLi too small");

    // Load
    let data_read: Vec<ReadNum> = file.dataset("data/read").expect("not found").read_raw().expect("E");
    let data_flag: Vec<u16> = file.dataset("data/flag").expect("not found").read_raw().expect("E");
    let data_rname_i: Vec<RNAMEi> = file.dataset("data/rname_i").expect("not found").read_raw().expect("E");
    let data_mapq: Vec<u8> = file.dataset("data/mapq").expect("not found").read_raw().expect("E");
    let data_cb_i: Vec<WLi> = file.dataset("data/cb_i").expect("not found").read_raw().expect("E");
    let data_ub_i: Vec<WLi> = file.dataset("data/ub_i").expect("not found").read_raw().expect("E");
    
    // Assert
    let mut set = HashSet::new(); for &value in &data_read { assert!(set.insert(value)); }
    assert_all_same(&[data_read.len(), data_flag.len(), data_rname_i.len(), data_mapq.len(), data_cb_i.len(), data_ub_i.len()]);
    
    // Zip
    let data = izip!(data_read, data_flag, data_rname_i, data_mapq, data_cb_i, data_ub_i)
               .map(|(i, f, r, m, c, u)| DataEntry{read: i, flag: f, rname_i: r, mapq: m, cb_i: c, ub_i: u})
               .collect::<Vec<DataEntry>>();
    
    return data;
}

/* SNV */
pub struct SNVEntry {
    pub read: ReadNum,
    pub pos: SeqPos,
    //pub ref: u8,
    pub alt: u8,
    pub qual: u8,
}
pub fn load_snv(file: &hdf5::File) -> Vec<SNVEntry> {
    assert!(file.dataset("snv/read").expect("not found").dtype().expect(".h5 error").size() <= size_of::<ReadNum>(), "ERROR: ReadNum too small");
    assert!(file.dataset("snv/pos").expect("not found").dtype().expect(".h5 error").size() <= size_of::<SeqPos>(), "ERROR: SeqPos too small");
    assert!(file.dataset("snv/alt").expect("not found").dtype().expect(".h5 error").size() <= size_of::<u8>(), "ERROR: u8 too small");
    assert!(file.dataset("snv/qual").expect("not found").dtype().expect(".h5 error").size() <= size_of::<u8>(), "ERROR: u8 too small");
    
    // Load
    let snv_read: Vec<ReadNum> = file.dataset("snv/read").expect("not found").read_raw().expect("E");
    let snv_pos: Vec<SeqPos> = file.dataset("snv/pos").expect("not found").read_raw().expect("E");
    let snv_alt: Vec<u8> = file.dataset("snv/alt").expect("not found").read_raw().expect("E");
    let snv_qual: Vec<u8> = file.dataset("snv/qual").expect("not found").read_raw().expect("E");
    
    // Assert
    assert_all_same(&[snv_read.len(), snv_pos.len(), snv_alt.len(), snv_qual.len()]);
    
    // Zip
    let snv = izip!(snv_read, snv_pos, snv_alt, snv_qual)
              .map(|(i, p, a, q)| SNVEntry{read: i, pos: p, alt: a, qual: q})
              .collect::<Vec<SNVEntry>>();
    
    return snv;
}

/* MATCHES */
pub struct MatchEntry {
    pub read: ReadNum,
    pub start: SeqPos,
    pub end: SeqPos,
}
pub fn load_matches(file: &hdf5::File) -> Vec<MatchEntry> {
    assert!(file.dataset("match/read").expect("not found").dtype().expect(".h5 error").size() <= size_of::<ReadNum>(), "ERROR: ReadNum too small");
    assert!(file.dataset("match/start").expect("not found").dtype().expect(".h5 error").size() <= size_of::<SeqPos>(), "ERROR: SeqPos too small");
    assert!(file.dataset("match/end").expect("not found").dtype().expect(".h5 error").size() <= size_of::<SeqPos>(), "ERROR: SeqPos too small");
    
    // Load
    let match_read: Vec<ReadNum> = file.dataset("match/read").expect("not found").read_raw().expect("E");
    let match_start: Vec<SeqPos> = file.dataset("match/start").expect("not found").read_raw().expect("E");
    let match_end: Vec<SeqPos> = file.dataset("match/end").expect("not found").read_raw().expect("E");

    // Assert
    assert_all_same(&[match_read.len(), match_start.len(), match_end.len()]);
    
    // Zip
    let matches = izip!(match_read, match_start, match_end)
                  .map(|(i, s, e)| MatchEntry{read: i, start: s, end: e})
                  .collect::<Vec<MatchEntry>>();
    
    return matches;
}

// Enables downsizing and writing to .h5
pub trait WriteVector {
    fn write_vector(&self, group: &hdf5::Group, name: &str);
}
impl<T: PrimInt> WriteVector for [T] {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        if self.len() == 0 {return ()}
        let max: u128 = self.iter().max().expect("ERROR: blank").to_u128().expect("ERROR: not u128");
        match max {
            0..=255            => group.new_dataset_builder().deflate(1).with_data(&self.into_iter().map(|x| x.to_u8().expect("ERROR: not u8")).collect::<Vec<u8>>()).create(name).expect(".h5 error"),
            256..=65535        => group.new_dataset_builder().deflate(1).with_data(&self.into_iter().map(|x| x.to_u16().expect("ERROR: not u16")).collect::<Vec<u16>>()).create(name).expect(".h5 error"),
            65536..=4294967295 => group.new_dataset_builder().deflate(1).with_data(&self.into_iter().map(|x| x.to_u32().expect("ERROR: not u32")).collect::<Vec<u32>>()).create(name).expect(".h5 error"),
            _                  => group.new_dataset_builder().deflate(1).with_data(&self.into_iter().map(|x| x.to_u64().expect("ERROR: not u64")).collect::<Vec<u64>>()).create(name).expect(".h5 error"),
        };
    }
}
