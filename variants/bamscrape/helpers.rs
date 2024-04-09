use std::collections::HashMap;
use num_traits::int::PrimInt;
use hdf5::types::VarLenAscii;
use csv::ReaderBuilder;
use std::fs::File;
use crate::WLi;

// create a (RNAME -> fasta byte offset) map from the .fai file
pub fn load_fai(index_path: &str) -> HashMap<String, usize> {
    let mut index_file = ReaderBuilder::new().delimiter(b'\t').has_headers(false)
                                             .from_reader(File::open(&index_path).expect("ERROR: could not open .fai file"));
    // loop through each record and save the fasta byte offset for each RNAME
    let mut map: HashMap<String, usize> = HashMap::new();
    for result in index_file.records() {
        let record = result.expect("ERROR: could not parse index file");
        if let (Some(key), Some(value)) = (record.get(0), record.get(2)) { // col 1 is RNAME, col 3 is fasta byte offset
            map.insert(key.to_string(), value.parse::<usize>().expect("ERROR: .fai column 3 could not be read as usize"));
        } else { panic!("ERROR: .fai file missing column data"); }
    }
    return map;
}

// helper method to assert that all elements of a vector are the same
pub fn assert_all_same(vec: &[usize]) {
    if let Some(first) = vec.first() {
        for element in vec.iter() {
            assert!(element == first, "ERROR: not all lengths are the same");
        }
    }
}

// Data structure that enables storing "indexes into a list of strings" instead of "strings"
pub struct Whitelist {
    map: HashMap<String, WLi>
}
impl Whitelist {
    // create an empty whitelist
    pub fn new() -> Self {
        Self { map: HashMap::new() }
    }
    // return the index of the string in the whitelist (add if needed)
    pub fn get(&mut self, string: &str) -> WLi {
        match self.map.get(string) {
            Some(&val) => val,
            None => {
                let n: WLi = self.map.len().try_into().expect("ERROR: WLi is too small");
                self.map.insert(string.to_string(), n);
                n
            }
        }
    }
    // return a vector of strings (in order)
    pub fn into_vector(self) -> Vec<String> {
        let mut pairs = self.map.into_iter().collect::<Vec<(String, WLi)>>();
        pairs.sort_by(|a, b| a.1.cmp(&b.1));
        pairs.into_iter().map(|(key, _)| key).collect::<Vec<String>>()
    }
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
impl WriteVector for Vec<String> {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        let vec: Vec<VarLenAscii> = self.into_iter().map(|s| VarLenAscii::from_ascii(s).expect("ERROR: non-ascii character")).collect();
        group.new_dataset_builder().with_data(&vec).create(name).expect(".h5 error");
    }
}
impl WriteVector for Vec<&str> {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        let vec: Vec<VarLenAscii> = self.into_iter().map(|s| VarLenAscii::from_ascii(s).expect("ERROR: non-ascii character")).collect();
        group.new_dataset_builder().with_data(&vec).create(name).expect(".h5 error");
    }
}

// assumptions:
// 2^64 is large enough to store any dataset
// xf only goes up to 63
// all stored values are UNSIGNED
