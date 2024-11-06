pub mod alphabet;
pub mod bwt;
pub mod fm_index;
pub mod occurrence;
pub mod search;
pub mod simd_instructions;
pub mod compressed_suffix_array;
//todo: how am I going to square the Suffix Array and the read collection?
//how to get the suffix array data out (use the sufr library its self)
//todo: count function, returns range, method on range to get len
//todo: set up hashmap for cannonical kmers
//todo: allow for adding 0-weighted references (linear walk through fastas)

fn create_sufr_suffix_array() -> Result<()> {}

#[cfg(test)]
mod tests {

    #[test]
    fn test_ceil() {
        assert_eq!((100 as usize).div_ceil(3), 34);
        assert_eq!((101 as usize).div_ceil(3), 34);
        assert_eq!((102 as usize).div_ceil(3), 34);
        assert_eq!((103 as usize).div_ceil(3), 35);
    }

    #[test]
    fn ascii_value_to_char(){
        assert_eq!(64 as char, 'A');
    }
}
