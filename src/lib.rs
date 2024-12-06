pub mod alphabet;
pub mod bwt;
pub mod compressed_suffix_array;
pub mod fm_index;
pub mod fm_index_file;
pub mod kmer_lookup_table;
pub mod search;
pub mod simd_instructions;
pub mod sequence_index;

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
    fn ascii_value_to_char() {
        assert_eq!(65 as char, 'A');
    }
}
