use serde::{Deserialize, Serialize};

/// Struct representing a sampled suffix array. Sampling the suffix array reduces the memory requirement, while still being able to reconstruct
/// the original position when used in conjunction with the rest of the FM-index.
/// The Suffix Array values stored inside this struct are bit-compressed to take up as little space as possible in memory.
///
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
pub(crate) struct CompressedSuffixArray {
    ///Actual Suffix Array values, compressed to remove leading zeros
    data: Vec<u64>,
    /// By what factor is the suffix array downsampled. a ratio of n only stores values whose indices are divisible by n.
    suffix_array_compression_ratio: usize,
    /// how many bits are required to store each value, i.e., an unsampled SA of length n will require log2_ceil(n) bits per element.
    bits_per_element: usize,
}

impl CompressedSuffixArray {
    /// Allocates space for a new Compressed Suffix Array, given the total uncompressed length and a compression ratio.
    ///
    pub(crate) fn new(uncompressed_length: usize, suffix_array_compression_ratio: usize) -> Self {
        assert_ne!(
            uncompressed_length, 0,
            "length of compressed suffix array should not be zero"
        );
        let largest_value_in_suffix_array = uncompressed_length - 1;
        let num_leading_zeros = largest_value_in_suffix_array.leading_zeros() as usize;
        let bits_per_element = 64 - num_leading_zeros;

        // find the ceil of the length of the suffix array divided by the compression ratio
        let num_compressed_elements = uncompressed_length.div_ceil(suffix_array_compression_ratio);
        let suffix_array_num_words =
            (num_compressed_elements as u128 * bits_per_element as u128).div_ceil(64) as usize;

        CompressedSuffixArray {
            data: vec![0; suffix_array_num_words],
            suffix_array_compression_ratio,
            bits_per_element,
        }
    }

    /// Returns a reference to the underlying compressed SA data
    ///
    pub(crate) fn data(&self) -> &Vec<u64> {
        &self.data
    }
    pub(crate) fn compression_ratio(&self) -> usize {
        self.suffix_array_compression_ratio
    }

    /// sets the value in the compressed suffix array.
    /// NOTE! the position is the compressed position, not the position in the full SA.
    ///
    pub(crate) fn set_value(&mut self, value: u64, position: usize) {
        let word_position = (position * self.bits_per_element as usize) / 64;
        let bit_position = (position * self.bits_per_element as usize) % 64;

        //create bitmasks to erase value here, and generate the values to write

        self.data[word_position] |= value << bit_position;
        if bit_position + self.bits_per_element > 64 {
            self.data[word_position + 1] |= match value.checked_shr(64 - bit_position as u32) {
                Some(val) => val,
                None => 0,
            };
        }
    }

    ///reconstructs the value at the given index in the suffix array.
    /// If that position wasn't sampled (i.e., not divisible by the compression ratio),
    /// this function will return None.
    pub(crate) fn reconstruct_value(&self, position: usize) -> Option<u64> {
        //if the position isn't sampled, return None to show it
        if position % self.suffix_array_compression_ratio as usize != 0 {
            return None;
        }
        let sampled_position = position / self.suffix_array_compression_ratio as usize;
        let word_position = (sampled_position * self.bits_per_element as usize) / 64;
        let first_word_bit_start_position =
            (sampled_position * self.bits_per_element as usize) % 64;
        let first_word_num_bits = self
            .bits_per_element
            .min(64 - first_word_bit_start_position);
        let second_word_num_bits = self.bits_per_element - first_word_num_bits;

        let first_word_bitmask = (1 << first_word_num_bits) - 1;
        let first_word_value_bits =
            (self.data[word_position] >> first_word_bit_start_position) & first_word_bitmask;
        //construct the second word in an IF statement because otherwise the array
        //could go out of bounds if the suffix array ends at this word.
        if second_word_num_bits != 0 {
            let second_word_bitmask = (1 << second_word_num_bits) - 1;
            let second_word_value_bits =
                (self.data[word_position + 1] & second_word_bitmask) << first_word_num_bits;

            return Some(first_word_value_bits | second_word_value_bits);
        } else {
            return Some(first_word_value_bits);
        }
    }

    /// Returns true if the given position is sampled in the compressed suffix array.
    pub fn position_is_sampled(&self, unsampled_position: usize) -> bool {
        return unsampled_position % self.suffix_array_compression_ratio == 0;
    }
}

#[cfg(test)]
mod tests {
    use super::CompressedSuffixArray;

    #[test]
    fn check_compressed_suffix_array() -> anyhow::Result<()> {
        let sa_len = 123451usize;
        let sa_values: Vec<u64> = (0..sa_len as u64).collect();

        for compression_ratio in 1usize..16 {
            let compressed_length = sa_len / compression_ratio;
            let mut csa = CompressedSuffixArray::new(sa_len as usize, compression_ratio);

            assert_eq!(
                compression_ratio,
                csa.compression_ratio(),
                "compression ratio did not match input"
            );

            //set the values
            for sa_value_idx in 0..compressed_length {
                csa.set_value(
                    sa_values[sa_value_idx * compression_ratio as usize],
                    sa_value_idx,
                );
            }

            //check the compressed suffix array values
            for sa_value_idx in 0..compressed_length {
                let expected_value = sa_values[sa_value_idx * compression_ratio as usize];
                let value_from_csa = csa
                    .reconstruct_value(sa_value_idx * compression_ratio as usize)
                    .expect("suffix array was not sampled here!");
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
