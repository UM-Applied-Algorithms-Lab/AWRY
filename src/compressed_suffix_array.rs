/// Struct representing a sampled suffix array. Sampling the suffix array reduces the memory requirement, while still being able to reconstruct
/// the original position when used in conjunction with the rest of the FM-index.
/// The Suffix Array values stored inside this struct are bit-compressed to take up as little space as possible in memory. 
pub struct CompressedSuffixArray {
    ///Actual Suffix Array values, compressed to remove leading zeros
    data: Vec<u64>, 
    ///Length of the uncompressed suffix array, this is different than data.len() 
    length: usize,
    /// By what factor is the suffix array downsampled. a ratio of n only stores values whose indices are divisible by n.
    suffix_array_compression_ratio: u64,
    /// how many bits are required to store each value, i.e., an unsampled SA of length n will require log2_ceil(n) bits per element. 
    bits_per_element: u64,
}

impl CompressedSuffixArray {
    /// Allocates space for a new Compressed Suffix Array, given the total uncompressed length and a compression ratio.
    pub fn new(length: usize, suffix_array_compression_ratio: u64) -> Self {
        //finds how many bits would be required to store any integer in the range 0 to length-1
        let mut num_bits_required: u64 = 0;
        let mut maximum_sa_value = length - 1;
        while maximum_sa_value != 0 {
            num_bits_required += 1;
            maximum_sa_value >>= 1;
        }

        CompressedSuffixArray {
            data: vec![0; length],
            length,
            suffix_array_compression_ratio,
            bits_per_element: num_bits_required,
        }
    }

    /// Returns a reference to the underlying compressed SA data
    pub fn data(&self) -> &Vec<u64> {
        &self.data
    }

    /// sets the value in the compressed suffix array.
    /// NOTE! the position is the compressed position, not the position in the full SA.
    pub fn set_value(&mut self, value: u64, position: usize) {
        let word_position = position / 64;
        let bit_position = position % 64;

        //create bitmasks to erase value here, and generate the values to write
        let first_word_bitmask = !(!0u64 << bit_position);
        let second_word_bitmask = !0u64 >> (64 - bit_position);
        let first_word_value = value << bit_position;
        let second_word_value = value >> (64 - bit_position);

        //apply the bitmasks and set the values
        self.data[word_position] &= first_word_bitmask;
        self.data[word_position + 1] &= second_word_bitmask;
        self.data[word_position] |= first_word_value;
        self.data[word_position + 1] |= second_word_value;
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
        let word_position = position / 64;
        let bit_position = position % 64;

        let unmasked_value =  self.data[word_position] >> bit_position
            | (self.data[word_position + 1] >> (64 - bit_position));
        let bitmask = (1<<self.bits_per_element)-1;

            return unmasked_value & bitmask;
    }

    /// Returns true if the given position is sampled in the compressed suffix array.
    pub fn position_is_sampled(&self, unsampled_position: u64) -> bool {
        return unsampled_position % self.suffix_array_compression_ratio == 0;
    }
}
