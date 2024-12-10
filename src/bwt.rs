use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    search::SearchPtr,
    simd_instructions::{Vec256, SimdVec256},
};
use serde::{Deserialize, Serialize};

pub const SIMD_ALIGNMENT_BYTES: usize = 32;

///block for a Nucleotide BWT. contains 6 milestones (packed to 8 for alignment), and 3 bit vectors
/// 
/// # Example
/// ```
/// use sufr_bwt::bwt::NucleotideBwtBlock;
/// 
/// let mut nucleotide_bwt_block = NucleotideBwtBlock::new();
/// nucleotide_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'), 0);
/// assert_eq!(nucleotide_bwt_block.ymbol_at(0), Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
/// ```
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
#[repr(align(32))]
pub struct NucleotideBwtBlock {
    bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
    milestones: [u64; Self::NUM_MILESTONES],
}

///block for a Amino BWT. contains 22 milestones (packed to 24 for alignment), and 5 bit vectors
/// 
/// # Example
/// ```
/// use sufr_bwt::bwt::AminoBwtBlock;
/// 
/// let mut amino_bwt_block = AminoBwtBlock::new();
/// amino_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'), 0);
/// assert_eq!(amino_bwt_block.symbol_at(0), Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'));
/// ``` 
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
#[repr(align(32))]
pub struct AminoBwtBlock {
    bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
    milestones: [u64; Self::NUM_MILESTONES],
}


impl NucleotideBwtBlock {
    pub const NUM_MILESTONES: usize = 8;
    pub const NUM_BIT_VECTORS: usize = 3;
    pub const NUM_SYMBOLS_PER_BLOCK: u64 = 256;

    ///creates a new bwt block, with data zeroed out.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// ``` 
    pub fn new() -> Self {
        NucleotideBwtBlock {
            bit_vectors: [Vec256::new(); Self::NUM_BIT_VECTORS],
            milestones: [0; Self::NUM_MILESTONES],
        }
    }

    ///creates a bwt block from the given milestone and bit-vector data.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::from_data(
    ///     [Vec256::new();NucleotideBwtBlock::NUM_BIT_VECTORS],
    ///     [1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64],
    /// );
    /// ```
    pub fn from_data(
        bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
        milestones: [u64; Self::NUM_MILESTONES],
    ) -> Self {
        NucleotideBwtBlock {
            bit_vectors,
            milestones,
        }
    }

    ///Returns the symbol at the given position in the BWT.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'), 0);
    /// assert_eq!(nucleotide_bwt_block.symbol_at(0), Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// ```
    pub fn symbol_at(&self, position_block: u64)->Symbol{
        let mut bit_vector_encoding: u64 = 0;

        for bit in 0..self.bit_vectors.len() {
            let bit_value = self.bit_vectors[bit].extract_bit(&position_block);
            bit_vector_encoding |= bit_value << bit;
        }

        Symbol::new_bit_vector(SymbolAlphabet::Nucleotide, bit_vector_encoding as u8)
    }

    ///Sets the symbol at the given position in the BWT block.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let mut nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'), 0);
    /// assert_eq!(nucleotide_bwt_block.symbol_at(0), Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// ```
    pub fn set_symbol_at(&mut self, symbol:&Symbol, position_in_block: u64){
        let mut encoded_symbol = symbol.bit_vector();

        //sets the bits in the bit-vectors based on the position and symbol given
        let mut bit_vector_idx = 0;
        while encoded_symbol != 0 {
            if encoded_symbol & 0x1 == 1 {
                self.bit_vectors[bit_vector_idx].set_bit(&position_in_block);
            }
            encoded_symbol >>= 1;
            bit_vector_idx += 1;
        }
    }

    ///sets this block's milestone values using the given vector
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let mut nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// assert_eq!(nucleotide_bwt_block.milestones(), &[1000, 2000, 3000, 4000]);
    /// ``` 
    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        debug_assert!(values.len() >= SymbolAlphabet::Nucleotide.cardinality() as usize);

        for milestone_idx in 0..SymbolAlphabet::Nucleotide.cardinality() as usize {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }

    ///gets this block's milestone corresponding to the given symbol.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// assert_eq!(nucleotide_bwt_block.milestone(&Symbol::new_index(SymbolAlphabet::Nucleotide, 0)), 1000);
    /// ```
    #[inline]
    pub fn milestone(&self, symbol: &Symbol) -> u64 {
        return self.milestones[symbol.index() as usize];
    }

    /// gets a reference to the milestones array
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// for milestone in nucleotide_bwt_block.milestones(){
    ///     println!("milestone: {}", milestone);
    /// }
    /// ```
    pub fn milestones(&self) -> &[u64] {
        &self.milestones
    }
    /// gets a reference to the bit_vectors array
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// for bit_vector in nucleotide_bwt_block.bit_vectors(){
    ///     println!("bit_vector: {:?}", bit_vector);
    /// }
    /// ```
    pub fn bit_vectors(&self) -> &[Vec256] {
        &self.bit_vectors
    }

    /// Gets the result of the occurrence function for the local position in this function.
    /// The occurrence function uses the milestone value and the masked occurrenc vector to
    /// determine how many instances of the given character were before this position.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::NucleotideBwtBlock;
    /// 
    /// let nucleotide_bwt_block = NucleotideBwtBlock::new();
    /// nucleotide_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// assert_eq!(nucleotide_bwt_block.global_occurrence(0, &Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A')), 0);
    /// ```
    #[inline]
    pub fn global_occurrence(&self, local_query_position: u64, symbol: &Symbol) -> u64 {
        let milestone_count = self.milestone(&symbol);
        let vec0 = SimdVec256::from(self.bit_vectors[0]);
        let vec1 = SimdVec256::from(self.bit_vectors[1]);
        let vec2 = SimdVec256::from(self.bit_vectors[2]);
        let occurrence_vector = match &symbol.index() {
            1 => vec2.and(&vec1), //A:    0b110
            2 => vec2.and(&vec0), //C:    0b101
            3 => vec1.and(&vec0), //G:    0b011
            4 => vec2.andnot(&vec0.andnot(&vec1)), //N:    0b010
            5 => vec2.andnot(&vec1.andnot(&vec0)), //T:    0b001
            _ => {
                panic!("illegal letter index given in global occurrence function symbol idx given: {}", symbol.index());
            } //assume every other character is an N, since it's illegal to search for a sentinel
        };

        let popcount = occurrence_vector.masked_popcount(local_query_position);

        return milestone_count + popcount as u64;
    }
}

impl AminoBwtBlock {
    pub const NUM_MILESTONES: usize = 24;
    pub const NUM_BIT_VECTORS: usize = 5;
    pub const NUM_SYMBOLS_PER_BLOCK: usize = 256;

    /// create a new bwt block, with data zeroed out
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::AminoBwtBlock;
    /// 
    /// let amino_bwt_block = AminoBwtBlock::new();
    /// ``` 
    pub fn new() -> Self {
        AminoBwtBlock {
            bit_vectors: [Vec256::new(); Self::NUM_BIT_VECTORS],
            milestones: [0; Self::NUM_MILESTONES],
        }
    }

    /// create a new bwt block from the given data.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::AminoBwtBlock;
    /// 
    /// let amino_bwt_block = AminoBwtBlock::from_data(
    ///     [Vec256::new();AminoBwtBlock::NUM_BIT_VECTORS],
    ///     [
    ///         1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64,
    ///         9000u64, 10000u64, 11000u64, 12000u64, 13000u64, 14000u64, 15000u64, 16000u64,
    ///         17000u64, 18000u64, 19000u64, 20000u64, 21000u64, 22000u64, 23000u64, 24000u64,
    ///     ],
    /// );
    /// ```
    pub fn from_data(
        bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
        milestones: [u64; Self::NUM_MILESTONES],
    ) -> Self {
        AminoBwtBlock {
            milestones,
            bit_vectors,
        }
    }

    ///Gets the symbol at the given position in the BWT block.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::AminoBwtBlock;
    /// 
    /// let amino_bwt_block = AminoBwtBlock::new();
    /// amino_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'), 0);
    /// assert_eq!(amino_bwt_block.symbol_at(0), Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'));    
    pub fn symbol_at(&self, position_block: u64)->Symbol{
        let mut bit_vector_encoding: u64 = 0;

        for bit in 0..self.bit_vectors.len() {
            let bit_value = self.bit_vectors[bit].extract_bit(&position_block);
            bit_vector_encoding |= bit_value << bit;
        }

        Symbol::new_bit_vector(SymbolAlphabet::Amino, bit_vector_encoding as u8)
    }

    ///Sets the symbol at the given position in the BWT block.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::AminoBwtBlock;
    /// 
    /// let mut amino_bwt_block = AminoBwtBlock::new();
    /// amino_bwt_block.set_symbol_at(&Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'), 0);
    /// assert_eq!(amino_bwt_block.symbol_at(0), Symbol::new_ascii(SymbolAlphabet::Amino, 'Q'));
    /// ``` 
    pub fn set_symbol_at(&mut self, symbol:&Symbol, position_in_block: u64){
        //create a bitmask, we'll use this to set the bit with an OR operation
        let mut encoded_symbol = symbol.bit_vector();

        //sets the bits in the bit-vectors based on the position and symbol given
        let mut bit_vector_idx = 0;
        while encoded_symbol != 0 {
            if encoded_symbol & 0x1 == 1 {
                self.bit_vectors[bit_vector_idx].set_bit(&position_in_block);
            }
            encoded_symbol >>= 1;
            bit_vector_idx += 1;
        }
    }


    /// sets the milestones for this block with the values given.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::AminoBwtBlock;
    /// 
    /// let mut amino_bwt_block = AminoBwtBlock::new();
    /// amino_bwt_block.set_milestones(&[1000, 2000, 3000, 4000]);
    /// assert_eq!(amino_bwt_block.milestones(), &[1000, 2000, 3000, 4000]);
    /// ```
    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        debug_assert!(values.len() >= SymbolAlphabet::Amino.cardinality() as usize);

        for milestone_idx in 0..SymbolAlphabet::Amino.cardinality() as usize {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }

    /// returns the milestone value corresponding to the given symbol
    #[inline]
    pub fn milestone(&self, symbol: &Symbol) -> u64 {
        return self.milestones[symbol.index() as usize];
    }

    /// returns a slice view of the milestones for this block
    pub fn milestones(&self) -> &[u64; Self::NUM_MILESTONES] {
        &self.milestones
    }

    /// returns a slice view of the bit_vectors for this block
    pub fn bit_vectors(&self) -> &[Vec256; Self::NUM_BIT_VECTORS] {
        &self.bit_vectors
    }

    /// Gets the result of the occurrence function for the local position in this function.
    /// The occurrence function uses the milestone value and the masked occurrenc vector to
    /// determine how many instances of the given character were before this position.
    #[inline]
    pub fn global_occurrence(&self, local_query_position: SearchPtr, symbol: &Symbol) -> SearchPtr {
        let milestone_count = self.milestone(symbol);
        let vec0 = SimdVec256::from(self.bit_vectors[0]);
        let vec1 = SimdVec256::from(self.bit_vectors[1]);
        let vec2 = SimdVec256::from(self.bit_vectors[2]);
        let vec3 = SimdVec256::from(self.bit_vectors[3]);   
        let vec4 = SimdVec256::from(self.bit_vectors[4]);
        let occurrence_vector = match symbol.index() {
            1 => vec3.and(&vec4.andnot(&vec2)), //A:    0b01100
            2 => vec3.andnot(&vec2).and(&vec1.and(&vec0)), //C:    0b10111
            3 => vec1.and(&vec4.andnot(&vec0)), //D:    0b00011
            4 => vec4.andnot(&vec2.and(&vec1)), //E: 0b00110
            5 => vec0.andnot(&vec3).and(&vec2.and(&vec1)), //F:    0b11110
            6 => vec2.andnot(&vec0.andnot(&vec4)), //G:    0b11010
            7 => vec2.andnot(&vec3).and(&vec1.and(&vec0)), //H: 0b11011
            8 => vec2.andnot(&vec1.andnot(&vec4)), //I:    0b11001
            9 => vec3.andnot(&vec1.andnot(&vec4)), //K:    0b10101
            10 => vec1.andnot(&vec0.andnot(&vec4)), //L:    0b11100
            11 => vec1.andnot(&vec3).and(&vec2.and(&vec0)), //M:    0b11101
            12 => vec0.or(&vec1).andnot(&vec2.andnot(&vec3)), //N:    0b01000
            13 => vec3.and(&vec4.andnot(&vec0)), //P:    0b01001,
            14 => vec3.or(&vec1).andnot(&vec0.andnot(&vec2)), //Q:    0b00100
            15 => vec3.andnot(&vec2.andnot(&vec4)), //R:    0b10011
            16 => vec3.and(&vec4.andnot(&vec1)), //S:    0b01010
            17 => vec2.and(&vec4.andnot(&vec0)), //T:    0b00101
            18 => vec3.andnot(&vec0.andnot(&vec4)), //V:    0b10110
            19 => vec3.or(&vec2).andnot(&vec1.andnot(&vec0)), //W:    0b00001
            20 => vec3.and(&vec2).and(&vec1.and(&vec0)), //Ambiguity character X:  0b11111
            21 => vec0.or(&vec2).andnot(&vec3.andnot(&vec1)), //Y:    0b00010
            // 0b00000 is sentinel, but since you can't search for sentinels, it is not included here.
            _ => {
                panic!("illegal letter index given in global occurrence function");
            }
        };

        let popcount = occurrence_vector.masked_popcount(local_query_position);

        return milestone_count + popcount as u64;
    }
}

/// enum representing a BWT, either as Nucleotide symbols or Amino symbols

#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
#[serde(untagged)]
pub enum Bwt {
    Nucleotide(Vec<NucleotideBwtBlock>),
    Amino(Vec<AminoBwtBlock>),
}



impl Bwt {
    pub const NUM_SYMBOLS_PER_BLOCK: u64 = 256;

    /// returns how many blocks are present in the bwt
    pub fn num_bwt_blocks(&self) -> usize {
        match self {
            Bwt::Nucleotide(vec) => vec.len(),
            Bwt::Amino(vec) => vec.len(),
        }
    }
    /// sets a single symbol at the given position in the bwt bit vectors.
    /// this function is meant to be run for every position in the bwt
    /// as a part of BWT data creation
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::{NucleotideBwtBlock, AminoBwtBlock, Bwt};
    /// 
    /// let mut bwt = Bwt::Nucleotide(vec![nucleotide_bwt_block]);
    /// bwt.set_symbol_at(&SearchPtr::new(0), &Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// assert_eq!(bwt.symbol_at(&SearchPtr::new(0)), Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// ``` 
    pub fn set_symbol_at(&mut self, bwt_position: &SearchPtr, symbol: &Symbol) {
        //find the block, byte, and bit of the data we're setting
        let bwt_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK;

        match self{
            Bwt::Nucleotide(vec) =>  vec[bwt_block_idx as usize].set_symbol_at(symbol, position_in_block),
            Bwt::Amino(vec) => vec[bwt_block_idx as usize].set_symbol_at(symbol, position_in_block),
        }
    }

    /// reconstructs the symbol stored at the given bwt position
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::{NucleotideBwtBlock, AminoBwtBlock, Bwt};
    /// 
    /// let mut bwt = Bwt::Nucleotide(vec![nucleotide_bwt_block]);
    /// bwt.set_symbol_at(&SearchPtr::new(0), &Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// assert_eq!(bwt.symbol_at(&SearchPtr::new(0)), Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'C'));
    /// ```
    pub fn symbol_at(&self, bwt_position: &SearchPtr) -> Symbol {
        //find the block, byte, and bit of the data we're setting
        let position_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK;

        let mut bit_vector_encoding: u64 = 0;
        match &self {
            Bwt::Nucleotide(vec) => {
                let alphabet = SymbolAlphabet::Nucleotide;
                let bwt_block = &vec[position_block_idx as usize];

                for bit in 0..bwt_block.bit_vectors.len() {
                    let bit_value = bwt_block.bit_vectors[bit].extract_bit(&position_in_block);
                    bit_vector_encoding |= bit_value << bit;
                }

                Symbol::new_bit_vector(alphabet, bit_vector_encoding as u8)
            }

            Bwt::Amino(vec) => {
                let alphabet = SymbolAlphabet::Amino;
                let bwt_block = &vec[position_block_idx as usize];

                for bit in 0..bwt_block.bit_vectors.len() {
                    let bit_value = bwt_block.bit_vectors[bit].extract_bit(&position_in_block);
                    bit_vector_encoding |= bit_value << bit;
                }

                Symbol::new_bit_vector(alphabet, bit_vector_encoding as u8)
            }
        }
    }

    ///sets the milestone values based on the given counts array
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::{NucleotideBwtBlock, AminoBwtBlock, Bwt};
    /// 
    /// let mut bwt = Bwt::Nucleotide(vec![nucleotide_bwt_block]);
    /// bwt.set_milestones(0, &[1000, 2000, 3000, 4000]);
    /// assert_eq!(bwt.milestone(&Symbol::new_index(SymbolAlphabet::Nucleotide, 0)), 1000);    
    /// ```
    pub fn set_milestones(&mut self, block_idx: usize, counts: &Vec<u64>) {
        match self {
            Bwt::Nucleotide(vec) => vec[block_idx].set_milestones(counts),
            Bwt::Amino(vec) => vec[block_idx].set_milestones(counts),
        }
    }

    /// finds the total occurrence value for the given symbol at the specified global position
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::bwt::{NucleotideBwtBlock, AminoBwtBlock, Bwt};
    /// 
    /// let mut bwt = Bwt::Nucleotide(vec![nucleotide_bwt_block]);
    /// bwt.set_milestones(0, &[1000, 2000, 3000, 4000]);
    /// assert_eq!(bwt.global_occurrence(0, &Symbol::new_index(SymbolAlphabet::Nucleotide, 0)), 0);    
    /// ```
    pub fn global_occurrence(
        &self,
        pointer_global_position: SearchPtr,
        symbol: &Symbol,
    ) -> SearchPtr {
        let block_idx: u64 = pointer_global_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let local_query_position: u64 = pointer_global_position % Self::NUM_SYMBOLS_PER_BLOCK;

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

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use rand::{Rng, SeedableRng};

    use crate::{alphabet::{Symbol, SymbolAlphabet}, simd_instructions::Vec256};
    use super::{AminoBwtBlock, NucleotideBwtBlock};

    #[test]
    fn mock_nucleotide_empty_bwt_block_test() {
        let mock_bwt_block = NucleotideBwtBlock::from_data(
            [Vec256::new();NucleotideBwtBlock::NUM_BIT_VECTORS],
            [
                1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64,
            ],
        );

        //make sure that with all zeros, each position and each symbol just returns the milestone
        for symbol_idx in 1..6 {
            for position in 0..256 {
                let occurrence_count = mock_bwt_block.global_occurrence(
                    position,
                    &Symbol::new_index(crate::alphabet::SymbolAlphabet::Nucleotide, symbol_idx),
                );

                assert_eq!(occurrence_count, mock_bwt_block.milestones[symbol_idx as usize], 
                    "nucleotide occurrence did not exactly match milestone in empty bwt block, count for sym {}, pos {}.", symbol_idx, position);
            }
        }
    }

    #[test]
    fn mock_nucleotide_preset_bwt_block_test(){
        //set up the comparison data
        let mut mock_bwt_block = NucleotideBwtBlock::from_data(
            [Vec256::new();NucleotideBwtBlock::NUM_BIT_VECTORS],
            [
                1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64,
            ],
            
        );
        let mut seeded_rng = rand::rngs::StdRng::seed_from_u64(2);

        let mut counts:HashMap<(u64,u8), u64> = HashMap::new();
        let mut current_counts:Vec<u64> = mock_bwt_block.milestones.to_vec();
        for position in 0..super::Bwt::NUM_SYMBOLS_PER_BLOCK{
            let symbol_idx = seeded_rng.gen_range(0..SymbolAlphabet::Nucleotide.cardinality()); 
            let symbol = &crate::alphabet::Symbol::new_index(SymbolAlphabet::Nucleotide, symbol_idx as u8);

            mock_bwt_block.set_symbol_at(symbol, position);

            //increment the current counts
            current_counts[symbol_idx as usize] += 1;

            //set the full counts table
            for idx in 0..SymbolAlphabet::Nucleotide.cardinality(){
                counts.insert((position, idx as u8), current_counts[idx as usize]);
            }
        }


        //check the global_occurrence results against the predicted results
        //make sure that with all zeros, each position and each symbol just returns the milestone
        for symbol_idx in 1..SymbolAlphabet::Nucleotide.cardinality() {
            for position in 0..256 {
                let occurrence_count = mock_bwt_block.global_occurrence(
                    position,
                    &Symbol::new_index(crate::alphabet::SymbolAlphabet::Nucleotide, symbol_idx),
                );
                let expected_value = match counts.get(&(position, symbol_idx)){
                    Some(val) => val,
                    None => {
                        println!("failed to get value from hash table at pos {} idx {}", position, symbol_idx);
                        panic!();
                    },
                };
                assert_eq!(occurrence_count, *expected_value, 
                    "nucleotide occurrence did not exactly match milestone in randomized bwt block, count for sym {}, pos {}.", symbol_idx, position);
            }
        }
    }

    #[test]
    fn mock_amino_empty_bwt_block_test() {
        let mock_bwt_block = AminoBwtBlock::from_data(
            [Vec256::new();AminoBwtBlock::NUM_BIT_VECTORS],
            [
                1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64,
                9000u64, 10000u64, 11000u64, 12000u64, 13000u64, 14000u64, 15000u64, 16000u64,
                17000u64, 18000u64, 19000u64, 20000u64, 21000u64, 22000u64, 23000u64, 24000u64,
            ],
        );

        //make sure that with all zeros, each position and each symbol just returns the milestone
        for symbol_idx in 1..6 {
            for position in 0..256 {
                let occurrence_count = mock_bwt_block.global_occurrence(
                    position,
                    &Symbol::new_index(crate::alphabet::SymbolAlphabet::Amino, symbol_idx),
                );

                assert_eq!(occurrence_count, mock_bwt_block.milestones[symbol_idx as usize], 
                    "amino occurrence did not exactly match milestone in empty bwt block, count for sym {}, pos {}.", symbol_idx, position);
            }
        }
    }

    #[test]
    fn mock_amino_preset_bwt_block_test(){
        //set up the comparison data
        let mut mock_bwt_block = AminoBwtBlock::from_data(
            [Vec256::new();AminoBwtBlock::NUM_BIT_VECTORS],
            [
                1000u64, 2000u64, 3000u64, 4000u64, 5000u64, 6000u64, 7000u64, 8000u64,
                9000u64, 10000u64, 11000u64, 12000u64, 13000u64, 14000u64, 15000u64, 16000u64,
                17000u64, 18000u64, 19000u64, 20000u64, 21000u64, 22000u64, 23000u64, 24000u64,
            ],
        );

        let mut counts:HashMap<(u64,u8), u64> = HashMap::new();
        let mut current_counts:Vec<u64> = mock_bwt_block.milestones.to_vec();
        let mut seeded_rng = rand::rngs::StdRng::seed_from_u64(6);
        for position in 0..super::Bwt::NUM_SYMBOLS_PER_BLOCK{
            let symbol_idx = seeded_rng.gen_range(0..SymbolAlphabet::Amino.cardinality()); 
            let symbol = &crate::alphabet::Symbol::new_index(SymbolAlphabet::Amino, symbol_idx as u8);

            mock_bwt_block.set_symbol_at(symbol, position);

            //increment the current counts
            current_counts[symbol_idx as usize] += 1;

            //set the full counts table
            for idx in 0..SymbolAlphabet::Amino.cardinality(){
                counts.insert((position, idx), current_counts[idx as usize]);
            }
        }


        //check the global_occurrence results against the predicted results
        //make sure that with all zeros, each position and each symbol just returns the milestone
        for symbol_idx in 1..SymbolAlphabet::Amino.cardinality() {
            for position in 0..256 {
                let occurrence_count = mock_bwt_block.global_occurrence(
                    position,
                    &Symbol::new_index(crate::alphabet::SymbolAlphabet::Amino, symbol_idx),
                );
                let expected_value = counts.get(&(position, symbol_idx)).expect("failed to get value from hash table");
                assert_eq!(occurrence_count, *expected_value, 
                    "amino occurrence did not exactly match milestone in randomized bwt block, count for sym {}, pos {}.", symbol_idx, position);
            }
        }
    }

}
