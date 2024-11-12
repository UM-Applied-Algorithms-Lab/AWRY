use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    search::SearchPtr,
    simd_instructions::SimdVec256,
};

#[derive(Clone)]
pub struct NucleotideBwtBlock {
    milestones: [u64; Self::NUM_MILESTONES],
    bit_vectors: [SimdVec256; Self::NUM_BIT_VECTORS],
}

#[derive(Clone)]
pub struct AminoBwtBlock {
    milestones: [u64; Self::NUM_MILESTONES],
    bit_vectors: [SimdVec256; Self::NUM_BIT_VECTORS],
}

impl NucleotideBwtBlock {
    pub const NUM_MILESTONES: usize = 8;
    pub const NUM_BIT_VECTORS: usize = 3;
    pub const NUM_SYMBOLS_PER_BLOCK: usize = 256;
    pub fn new() -> Self {
        NucleotideBwtBlock {
            milestones: [0; Self::NUM_MILESTONES],
            bit_vectors: [SimdVec256::zero(); Self::NUM_BIT_VECTORS],
        }
    }

    pub fn from_data(
        milestones: [u64; Self::NUM_MILESTONES],
        bit_vectors: [SimdVec256; Self::NUM_BIT_VECTORS],
    ) -> Self {
        NucleotideBwtBlock {
            milestones,
            bit_vectors,
        }
    }

    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        for milestone_idx in 0..Self::NUM_MILESTONES {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }

    #[inline]
    pub fn get_milestone(&self, symbol: &Symbol) -> u64 {
        return self.milestones[symbol.index() as usize];
    }
    pub fn milestones(&self) -> [u64; Self::NUM_MILESTONES] {
        self.milestones
    }
    pub fn bit_vectors(&self) -> [SimdVec256; Self::NUM_BIT_VECTORS] {
        self.bit_vectors
    }

    #[inline]
    pub fn global_occurrence(&self, local_query_position: u64, symbol: &Symbol) -> u64 {
        let milestone_count = self.get_milestone(&symbol);
        let vectors = &self.bit_vectors;
        let occurrence_vector = match &symbol.index() {
            0 => vectors[2].and(&vectors[1]), //A:    0b110
            1 => vectors[2].and(&vectors[0]), //C:    0b101
            2 => vectors[1].and(&vectors[0]), //G:    0b011
            3 => vectors[2].andnot(&vectors[0].andnot(&vectors[1])), //N:    0b010
            4 => vectors[2].andnot(&vectors[1].andnot(&vectors[0])), //T:    0b001
            _ => {
                panic!("illegal letter index given in global occurrence function");
            } //assume every other character is an N, since it's illegal to search for a sentinel
        };

        let popcount = occurrence_vector.masked_popcount(local_query_position);

        //generate the position bitmask

        return milestone_count + popcount as u64;
    }
}

impl AminoBwtBlock {
    pub const NUM_MILESTONES: usize = 24;
    pub const NUM_BIT_VECTORS: usize = 5;
    pub const NUM_SYMBOLS_PER_BLOCK: usize = 256;
    pub fn new() -> Self {
        AminoBwtBlock {
            milestones: [0; Self::NUM_MILESTONES],
            bit_vectors: [SimdVec256::zero(); Self::NUM_BIT_VECTORS],
        }
    }
    pub fn from_data(
        milestones: [u64; Self::NUM_MILESTONES],
        bit_vectors: [SimdVec256; Self::NUM_BIT_VECTORS],
    ) -> Self {
        AminoBwtBlock {
            milestones,
            bit_vectors,
        }
    }
    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        for milestone_idx in 0..Self::NUM_MILESTONES {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }
    #[inline]
    pub fn get_milestone(&self, symbol: &Symbol) -> u64 {
        return self.milestones[symbol.index() as usize];
    }

    pub fn milestones(&self) -> [u64; Self::NUM_MILESTONES] {
        self.milestones
    }

    pub fn bit_vectors(&self) -> [SimdVec256; Self::NUM_BIT_VECTORS] {
        self.bit_vectors
    }
    #[inline]
    pub fn global_occurrence(&self, local_query_position: SearchPtr, symbol: &Symbol) -> SearchPtr {
        let milestone_count = self.get_milestone(symbol);
        let vecs = &self.bit_vectors;
        let occurrence_vector = match symbol.index() {
            0 => vecs[3].and(&vecs[4].andnot(&vecs[2])), //A:    0b01100
            1 => vecs[3].andnot(&vecs[2]).and(&vecs[1].andnot(&vecs[0])), //C:    0b10111
            2 => vecs[1].and(&vecs[4].andnot(&vecs[0])), //D:    0b00011
            3 => vecs[4].andnot(&vecs[2].and(&vecs[1])), //E: 0b00110
            4 => vecs[0].andnot(&vecs[3]).and(&vecs[2].and(&vecs[1])), //F:    0b11110
            5 => vecs[2].andnot(&vecs[0].andnot(&vecs[4])), //G:    0b11010
            6 => vecs[2].andnot(&vecs[3]).and(&vecs[1].and(&vecs[0])), //H: 0b11011
            7 => vecs[2].andnot(&vecs[1].andnot(&vecs[4])), //I:    0b11001
            8 => vecs[3].andnot(&vecs[1].andnot(&vecs[4])), //K:    0b10101
            9 => vecs[1].andnot(&vecs[0].andnot(&vecs[4])), //L:    0b11100
            10 => vecs[1].andnot(&vecs[3]).and(&vecs[2].and(&vecs[0])), //M:    0b11101
            11 => vecs[0].or(&vecs[1]).andnot(&vecs[2].andnot(&vecs[3])), //N:    0b01000
            12 => vecs[3].and(&vecs[4].andnot(&vecs[0])), //P:    0b01001,
            13 => vecs[3].or(&vecs[1]).andnot(&vecs[0].andnot(&vecs[2])), //Q:    0b00100
            14 => vecs[3].andnot(&vecs[2].andnot(&vecs[4])), //R:    0b10011
            15 => vecs[3].and(&vecs[4].andnot(&vecs[1])), //S:    0b01010
            16 => vecs[2].and(&vecs[4].andnot(&vecs[0])), //T:    0b00101
            17 => vecs[3].andnot(&vecs[0].andnot(&vecs[4])), //V:    0b10110
            18 => vecs[3].or(&vecs[2]).andnot(&vecs[1].andnot(&vecs[0])), //W:    0b00001
            19 => vecs[3].and(&vecs[2]).and(&vecs[1].and(&vecs[0])), //Ambiguity character X:  0b11111
            20 => vecs[0].or(&vecs[2]).andnot(&vecs[3].andnot(&vecs[1])), //Y:    0b00010
            // 0b00000 is sentinel, but since you can't search for sentinels, it is not included here.
            _ => {
                panic!("illegal letter index given in global occurrence function");
            }
        };

        let popcount = occurrence_vector.masked_popcount(local_query_position);
        return milestone_count + popcount as u64;
    }
}

pub enum Bwt {
    Nucleotide(Vec<NucleotideBwtBlock>),
    Amino(Vec<AminoBwtBlock>),
}

impl Bwt {
    pub const NUM_SYMBOLS_PER_BLOCK: usize = 256;
    pub fn num_bwt_blocks(&self) -> usize {
        match self {
            Bwt::Nucleotide(vec) => vec.len(),
            Bwt::Amino(vec) => vec.len(),
        }
    }
    /// sets a single symbol at the given position in the bwt bit vectors.
    /// this function is meant to be run for every position in the bwt
    /// as a part of BWT data creation
    pub fn set_symbol_at(&self, bwt_position: &SearchPtr, symbol: &Symbol) {
        //find the block, byte, and bit of the data we're setting
        let position_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK as u64;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK as u64;

        //create a bitmask, we'll use this to set the bit with an OR operation
        let vector_bitmask = SimdVec256::as_one_hot(position_in_block);
        let mut encoded_symbol = symbol.bit_vector();

        let mut bit_vector_idx = 0;
        while encoded_symbol != 0 {
            if encoded_symbol & 0x1 == 0 {
                match self {
                    Bwt::Nucleotide(vec) => {
                        vec[position_block_idx as usize].bit_vectors[bit_vector_idx]
                            .or(&vector_bitmask);
                    }
                    Bwt::Amino(vec) => {
                        vec[position_block_idx as usize].bit_vectors[bit_vector_idx]
                            .or(&vector_bitmask);
                    }
                }
            }
            encoded_symbol >>= 1;
            bit_vector_idx += 1;
        }
    }

    pub fn get_symbol_at(&self, bwt_position: &SearchPtr) -> Symbol {
        //find the block, byte, and bit of the data we're setting
        let position_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK as u64;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK as u64;

        let mut bit_vector_encoding: u64 = 0;
        match &self {
            Bwt::Nucleotide(vec) => {
                let alphabet = SymbolAlphabet::Nucleotide;
                let bwt_block = &vec[position_block_idx as usize];

                for bit in 0..bwt_block.bit_vectors.len() {
                    let bit_value = bwt_block.bit_vectors[bit].get_bit(&position_in_block);
                    bit_vector_encoding |= bit_value << bit;
                }

                Symbol::new_bit_vector(alphabet, bit_vector_encoding as u8)
            }

            Bwt::Amino(vec) => {
                let alphabet = SymbolAlphabet::Nucleotide;
                let bwt_block = &vec[position_block_idx as usize];

                for bit in 0..bwt_block.bit_vectors.len() {
                    let bit_value = bwt_block.bit_vectors[bit].get_bit(&position_in_block);
                    bit_vector_encoding |= bit_value << bit;
                }

                Symbol::new_bit_vector(alphabet, bit_vector_encoding as u8)
            }
        }
    }

    pub fn set_milestones(&mut self, block_idx: usize, counts: &Vec<u64>) {
        match self {
            Bwt::Nucleotide(vec) => vec[block_idx].set_milestones(counts),
            Bwt::Amino(vec) => vec[block_idx].set_milestones(counts),
        }
    }
    fn get_milestone(&self, block_idx: usize, letter_idx: u8) -> u64 {
        return match self {
            Bwt::Nucleotide(vec) => vec[block_idx].milestones[letter_idx as usize],
            Bwt::Amino(vec) => vec[block_idx].milestones[letter_idx as usize],
        };
    }
    pub fn global_occurrence(
        &self,
        pointer_global_position: SearchPtr,
        symbol: &Symbol,
    ) -> SearchPtr {
        let block_idx: u64 = pointer_global_position / Self::NUM_SYMBOLS_PER_BLOCK as u64;
        let local_query_position: u64 =
            pointer_global_position % Self::NUM_SYMBOLS_PER_BLOCK as u64;

        match self {
            Bwt::Nucleotide(vec) => {
                vec[block_idx as usize].global_occurrence(local_query_position, symbol)
            }

            Bwt::Amino(vec) => {
                vec[block_idx as usize].global_occurrence(local_query_position, symbol)
            }
        }
    }
}
