use std::collections::HashMap;
use num_traits::int::PrimInt;
use hdf5::types::VarLenAscii;
use crate::WLi;

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
                self.map.insert(string.to_owned(), n);
                n
            }
        }
    }
    // return a vector of strings (in order)
    pub fn into_vector(self) -> Vec<String> {
        let mut pairs = self.map.into_iter().collect::<Vec<(String, WLi)>>();
        pairs.sort_by(|a, b| a.1.cmp(&b.1));
        assert!(pairs.iter().enumerate().all(|(i1, &(_, i2))| i1 == i2 as usize));
        pairs.into_iter().map(|(key, _)| key).collect::<Vec<String>>()
    }
}

// Enables downsizing and writing to .h5
pub trait WriteVector {
    fn write_vector(&self, group: &hdf5::Group, name: &str);
}
impl<T: PrimInt> WriteVector for [T] {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        let max: u128 = self.iter().max().expect("ERROR: blank").to_u128().expect("ERROR: not u128");
        match max {
            0..=255            => group.new_dataset_builder().deflate(1).with_data(&self.iter().map(|x| x.to_u8().expect("ERROR: not u8")).collect::<Vec<u8>>()).create(name).expect(".h5 error"),
            256..=65535        => group.new_dataset_builder().deflate(1).with_data(&self.iter().map(|x| x.to_u16().expect("ERROR: not u16")).collect::<Vec<u16>>()).create(name).expect(".h5 error"),
            65536..=4294967295 => group.new_dataset_builder().deflate(1).with_data(&self.iter().map(|x| x.to_u32().expect("ERROR: not u32")).collect::<Vec<u32>>()).create(name).expect(".h5 error"),
            _                  => group.new_dataset_builder().deflate(1).with_data(&self.iter().map(|x| x.to_u64().expect("ERROR: not u64")).collect::<Vec<u64>>()).create(name).expect(".h5 error"),
        };
    }
}
impl<T: AsRef<[u8]>> WriteVector for Vec<T> {
    fn write_vector(&self, group: &hdf5::Group, name: &str) {
        let vec: Vec<VarLenAscii> = self.iter().map(|s| VarLenAscii::from_ascii(s).expect("ERROR: non-ascii character")).collect();
        group.new_dataset_builder().deflate(1).with_data(&vec).create(name).expect(".h5 error");
    }
}
