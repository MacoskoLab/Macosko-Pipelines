use std::collections::HashMap;
use std::collections::HashSet;
use num_traits::int::PrimInt;
use itertools::izip;

// helper method to assert that all elements of a vector are unique
fn assert_unique(values: &[u64]) {
    let mut set = HashSet::new();
    for &value in values {
        assert!(set.insert(value));
    }
}
// helper method to get the average value of a list of u8 (used for computing average mapq)
pub fn average_round_u8(values: &[u8]) -> u8 {
    let sum: u32 = values.iter().map(|&val| u32::from(val)).sum();
    let len = values.len() as u32;
    let avg = (sum + len/2) / len;
    return avg as u8
}
// helper method to assert that all elements of a vector are the same
pub fn assert_all_same(vec: &[usize]) {
    if let Some(first) = vec.first() {
        for element in vec.iter() {
            assert!(element == first, "ERROR: not all lengths are the same");
        }
    }
}

pub type ReadNum = u32;
pub type SeqPos = u32;
pub type WLi = u32;
pub type RNAMEi = u16;

#[derive(Debug)]
pub struct DataEntry {
    pub read: ReadNum,
    pub rname_i: RNAMEi,
    pub mapq: u8,
    pub cb_i: WLi,
    pub ub_i: WLi,
}
pub fn load_data(file: &hdf5::File) -> Vec<DataEntry> {
    // Load
    let data_read: Vec<u64> = file.dataset("data/read").expect("not found").read_raw().expect("E"); assert_unique(&data_read);
    let data_rname_i: Vec<u64> = file.dataset("data/rname_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_rname_i.len());
    let data_mapq: Vec<u8> = file.dataset("data/mapq").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_mapq.len());
    let data_cb_i: Vec<u64> = file.dataset("data/cb_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_cb_i.len());
    let data_ub_i: Vec<u64> = file.dataset("data/ub_i").expect("not found").read_raw().expect("E"); assert!(data_read.len() == data_ub_i.len());
    
    // Downsize
    assert!(data_read.iter().max().expect("E") <= &(ReadNum::MAX as u64));
    assert!(data_rname_i.iter().max().expect("E") <= &(RNAMEi::MAX as u64));
    // assert!(data_mapq.iter().max().expect("E") <= &(u8::MAX as u64));
    assert!(data_cb_i.iter().max().expect("E") <= &(WLi::MAX as u64));
    assert!(data_ub_i.iter().max().expect("E") <= &(WLi::MAX as u64));
    let mut data = izip!(data_read, data_rname_i, data_mapq, data_cb_i, data_ub_i)
    .map(|(a, b, c, d, e)| DataEntry{read: a as ReadNum, rname_i: b as RNAMEi, mapq: c as u8, cb_i: d as WLi, ub_i: e as WLi})
    .collect::<Vec<DataEntry>>();
    
    // Sort by cb, umi
    data.sort_by(|d1, d2| {
        (d1.cb_i).cmp(&d2.cb_i).then((d1.ub_i).cmp(&d2.ub_i))
    });
    
    return data;
}

#[derive(Debug)]
pub struct SNVEntry {
    pub read: ReadNum,
    pub pos: SeqPos,
    pub alt: u8,
    pub qual: u8,
}
pub fn load_snv(file: &hdf5::File, order_map: &HashMap<ReadNum, usize>) -> Vec<SNVEntry> {
    // Load
    let snv_read: Vec<u64> = file.dataset("snv/read").expect("not found").read_raw().expect("E");
    let snv_pos: Vec<u64> = file.dataset("snv/pos").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_pos.len());
    let snv_alt: Vec<u8> = file.dataset("snv/alt").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_alt.len());
    let snv_qual: Vec<u8> = file.dataset("snv/qual").expect("not found").read_raw().expect("E"); assert!(snv_read.len() == snv_qual.len());
    
    // Downsize
    assert!(snv_read.iter().max().expect("E") <= &(ReadNum::MAX as u64));
    assert!(snv_pos.iter().max().expect("E") <= &(SeqPos::MAX as u64));
    // assert!(snv_alt.iter().max().expect("E") <= &(u8::MAX as u64));
    // assert!(snv_qual.iter().max().expect("E") <= &(u8::MAX as u64));
    let mut snv = izip!(snv_read, snv_pos, snv_alt, snv_qual)
    .map(|(a, b, c, d)| SNVEntry{read: a as ReadNum, pos: b as SeqPos, alt: c as u8, qual: d as u8})
    .collect::<Vec<SNVEntry>>();
    
    // Sort by umi
    snv.sort_by_key(|entry| order_map.get(&entry.read).expect("E"));
    
    return snv;
}

#[derive(Debug)]
pub struct MatchEntry {
    pub read: ReadNum,
    pub start: SeqPos,
    pub end: SeqPos,
}
pub fn load_matches(file: &hdf5::File, order_map: &HashMap<ReadNum, usize>) -> Vec<MatchEntry> {
    // Load
    let match_read: Vec<u64> = file.dataset("match/read").expect("not found").read_raw().expect("E");
    let match_start: Vec<u64> = file.dataset("match/start").expect("not found").read_raw().expect("E"); assert!(match_read.len() == match_start.len());
    let match_end: Vec<u64> = file.dataset("match/end").expect("not found").read_raw().expect("E"); assert!(match_read.len() == match_end.len());

    // Downsize
    assert!(match_read.iter().max().expect("E") <= &(ReadNum::MAX as u64));
    assert!(match_start.iter().max().expect("E") <= &(SeqPos::MAX as u64));
    assert!(match_end.iter().max().expect("E") <= &(SeqPos::MAX as u64));
    let mut matches = izip!(match_read, match_start, match_end)
    .map(|(a, b, c)| MatchEntry{read: a as ReadNum, start: b as SeqPos, end: c as SeqPos})
    .collect::<Vec<MatchEntry>>();
    
    // Sort by umi
    matches.sort_by_key(|entry| order_map.get(&entry.read).expect("E"));
    
    return matches;
}

// Data structure that enables storing "indexes into a list of SNV" instead of "SNV"
pub struct SNVTable {
    map: HashMap<(RNAMEi,SeqPos,u8), usize>
}
impl SNVTable {
    // create an empty whitelist
    pub fn new() -> Self {
        Self { map: HashMap::new() }
    }
    // return the index of the string in the whitelist (add if needed)
    pub fn get(&mut self, snv: &(RNAMEi,SeqPos,u8)) -> usize {
        match self.map.get(snv) {
            Some(&val) => val,
            None => {
                let n: usize = self.map.len();
                self.map.insert(*snv, n);
                n
            }
        }
    }
    // return a vector of tuples (in order)
    pub fn into_vectors(self) -> (Vec<usize>,Vec<RNAMEi>,Vec<SeqPos>,Vec<u8>) {
        let mut pairs: Vec<((RNAMEi,SeqPos,u8), usize)> = self.map.into_iter().collect();
        pairs.sort_by(|a, b| a.1.cmp(&b.1));
        let mut vec_snv_i: Vec<usize> = Vec::new();
        let mut vec_rname_i: Vec<RNAMEi> = Vec::new();
        let mut vec_pos: Vec<SeqPos> = Vec::new();
        let mut vec_alt: Vec<u8> = Vec::new();
        for ((rname_i, pos, alt), snv_i) in pairs {
            vec_snv_i.push(snv_i);
            vec_rname_i.push(rname_i);
            vec_pos.push(pos);
            vec_alt.push(alt);
        }
        (vec_snv_i,vec_rname_i, vec_pos, vec_alt)
    }
}

// Data structure that aggregates all the information of a UMI
#[derive(Debug)]
pub struct UMIinfo {
    pub reads: ReadNum,
    pub mapq: Vec<u8>,
    pub snv_i: HashMap<(usize,SeqPos),(ReadNum,ReadNum)>, // (snv_i, pos) -> (hq, lq)
    pub start: Vec<SeqPos>,
    pub end: Vec<SeqPos>,
}
impl UMIinfo {
    pub fn new() -> Self {
        Self { reads: 0, mapq: Vec::new(), snv_i: HashMap::new(), start: Vec::new(), end: Vec::new(), }
    }
}

// given a series of ranges, collapse into a consensus range
pub fn merge_intervals(starts: &Vec<SeqPos>, ends: &Vec<SeqPos>) -> Vec<(SeqPos,SeqPos)> {
    let mut intervals: Vec<(SeqPos,SeqPos)> = starts.iter().zip(ends.iter()).map(|(&start, &end)| (start, end)).collect();
    intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut merged_intervals: Vec<(SeqPos,SeqPos)> = Vec::new();
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

// Enables downsizing and writing to .h5
pub trait WriteVector {
    fn write_vector(&self, group: &hdf5::Group, name: &str);
}
impl<T: PrimInt> WriteVector for [T] {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        let max: u64 = self.iter().max().expect("ERROR: blank").to_u64().expect("ERROR: not u64");
        match max {
            0..=255            => group.new_dataset_builder().with_data(&self.into_iter().map(|x| x.to_u8().expect("ERROR: not u8")).collect::<Vec<u8>>()).create(name).expect(".h5 error"),
            256..=65535        => group.new_dataset_builder().with_data(&self.into_iter().map(|x| x.to_u16().expect("ERROR: not u16")).collect::<Vec<u16>>()).create(name).expect(".h5 error"),
            65536..=4294967295 => group.new_dataset_builder().with_data(&self.into_iter().map(|x| x.to_u32().expect("ERROR: not u32")).collect::<Vec<u32>>()).create(name).expect(".h5 error"),
            _                  => group.new_dataset_builder().with_data(&self.into_iter().map(|x| x.to_u64().expect("ERROR: not u64")).collect::<Vec<u64>>()).create(name).expect(".h5 error"),
        };
    }
}
