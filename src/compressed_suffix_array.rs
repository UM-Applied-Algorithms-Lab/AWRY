/// Struct representing a sampled suffix array. Sampling the suffix array reduces the memory requirement, while still being able to reconstruct
/// the original position when used in conjunction with the rest of the FM-index.
/// The Suffix Array values stored inside this struct are bit-compressed to take up as little space as possible in memory.
pub struct CompressedSuffixArray {
    ///Actual Suffix Array values, compressed to remove leading zeros
    data: Vec<u64>,
    /// By what factor is the suffix array downsampled. a ratio of n only stores values whose indices are divisible by n.
    suffix_array_compression_ratio: u64,
    /// how many bits are required to store each value, i.e., an unsampled SA of length n will require log2_ceil(n) bits per element.
    bits_per_element: u64,
}

impl CompressedSuffixArray {
    /// Allocates space for a new Compressed Suffix Array, given the total uncompressed length and a compression ratio.
    pub fn new(length: usize, suffix_array_compression_ratio: u64) -> Self {
        let bits_per_element = 64 - (length - 1).leading_zeros() as u64;

        CompressedSuffixArray {
            data: vec![0; length],
            suffix_array_compression_ratio,
            bits_per_element,
        }
    }

    /// Returns a reference to the underlying compressed SA data
    pub fn data(&self) -> &Vec<u64> {
        &self.data
    }
    pub fn compression_ratio(&self) -> u64 {
        self.suffix_array_compression_ratio
    }

    /// sets the value in the compressed suffix array.
    /// NOTE! the position is the compressed position, not the position in the full SA.
    pub fn set_value(&mut self, value: u64, position: usize) {
        let word_position = (position * self.bits_per_element as usize) / 64;
        let bit_position = (position * self.bits_per_element as usize) % 64;

        //create bitmasks to erase value here, and generate the values to write

        self.data[word_position] |= value << bit_position;
        self.data[word_position + 1] |= match value.checked_shr(64 - bit_position as u32) {
            Some(val) => val,
            None => 0,
        };
    }

    pub fn get_value(&self, position: usize) -> Option<u64> {
        if position % self.suffix_array_compression_ratio as usize == 0 {
            Some(self.reconstruct_value(position / self.suffix_array_compression_ratio as usize))
        } else {
            None
        }
    }

    ///reconstructs the value stored in the given element position in the compressed suffix array.
    /// Note that this position does not represent the uncompressed SA position, but instead the
    /// ith index in the compressed SA.
    fn reconstruct_value(&self, position: usize) -> u64 {
        let word_position = (position * self.bits_per_element as usize) / 64;
        let bit_position = (position * self.bits_per_element as usize) % 64;
        let bitmask = (1 << self.bits_per_element) - 1;
        let first_word = (self.data[word_position] & (bitmask << bit_position)) >> bit_position;
        let second_word = match self.data[word_position + 1].checked_shr(64 - bit_position as u32) {
            Some(val) => val & (bitmask >> 64 - bit_position),
            None => 0,
        };
        let unmasked_value = first_word | second_word;

        return unmasked_value;
    }

    /// Returns true if the given position is sampled in the compressed suffix array.
    pub fn position_is_sampled(&self, unsampled_position: u64) -> bool {
        return unsampled_position % self.suffix_array_compression_ratio == 0;
    }
}

#[cfg(test)]
mod tests {
    use super::CompressedSuffixArray;

    #[test]
    fn check_compressed_suffix_array() -> anyhow::Result<()> {
        let sa_len = 256;
        let sa_values: Vec<u64> = (0..sa_len).collect();
        
        for compression_ratio in 1..16 {
            let compressed_length = (sa_len / compression_ratio) as usize;
            let mut csa = CompressedSuffixArray::new(sa_len as usize, compression_ratio);

            assert_eq!(
                compression_ratio,
                csa.compression_ratio(),
                "compression ratio did not match input"
            );

            //set the values
            for sa_value_idx in 0..compressed_length{
                csa.set_value(
                    sa_values[sa_value_idx * compression_ratio as usize],
                    sa_value_idx,
                );
            }

            //check the compressed suffix array values
            for sa_value_idx in 0..compressed_length{
                let expected_value = sa_values[sa_value_idx * compression_ratio as usize];
                let value_from_csa = csa.reconstruct_value(sa_value_idx);
                assert_eq!(
                    value_from_csa,
                    expected_value,
                    "sa_idx {} failed for csr {}, supposed to be {}, found {}. bits per is {}",
                    sa_value_idx,
                    compression_ratio,
                    expected_value,
                    value_from_csa,
                    csa.bits_per_element
                );
            }
        }

        Ok(())
    }

    #[test]
    fn check_bits_per_element() -> anyhow::Result<()> {
        for (length, expected_bits) in [
            (15, 4),
            (16, 4),
            (17, 5),
            (31, 5),
            (32, 5),
            (33, 6),
            (1022, 10),
            (1023, 10),
            (1024, 10),
            (1025, 11),
            (65535, 16),
            (65536, 16),
            (65537, 17),
            (2usize.pow(31) - 1, 31),
            (2usize.pow(31), 31),
            (2usize.pow(31) + 1, 32),
        ] {
            let csa = CompressedSuffixArray::new(length, 8);

            assert_eq!(
                expected_bits, csa.bits_per_element,
                "bits did not match expected for length {}, bits {}",
                length, expected_bits
            );
        }

        Ok(())
    }
}
