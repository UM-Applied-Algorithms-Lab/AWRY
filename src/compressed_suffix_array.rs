
pub struct CompressedSuffixArray {
    data: Vec<u64>,
    length: usize,
    suffix_array_compression_ratio: usize,
    bits_per_element: u64,
}

impl CompressedSuffixArray {
    pub fn new(length: usize, suffix_array_compression_ratio: usize) -> Self {
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
        if position % self.suffix_array_compression_ratio == 0 {
            Some(self.reconstruct_value(position / self.suffix_array_compression_ratio))
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

        return self.data[word_position] >> bit_position
            | (self.data[word_position + 1] >> (64 - bit_position));
    }
}
