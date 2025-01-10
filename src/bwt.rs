

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    search::SearchPtr,
    simd_instructions::{masked_popcount, simd_and, simd_andnot, simd_or, Vec256},
};
use mem_dbg::MemSize;
use serde::{Deserialize, Serialize};

///block for a Nucleotide BWT. contains 6 milestones (packed to 8 for alignment), and 3 bit vectors
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default, MemSize)]
#[repr(align(32))]
pub (crate) struct NucleotideBwtBlock {
    bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
    milestones: [u64; Self::NUM_MILESTONES],
}

///block for a Amino BWT. contains 22 milestones (packed to 24 for alignment), and 5 bit vectors
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default, MemSize)]
#[repr(align(32))]
pub  (crate) struct AminoBwtBlock {
    bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
    milestones: [u64; Self::NUM_MILESTONES],
}


impl NucleotideBwtBlock {
    pub (crate) const NUM_MILESTONES: usize = 8;
    pub (crate) const NUM_BIT_VECTORS: usize = 3;

    ///creates a new bwt block, with data zeroed out.

    pub (crate) fn new() -> Self {
        NucleotideBwtBlock {
            bit_vectors: [Vec256::new(); Self::NUM_BIT_VECTORS],
            milestones: [0; Self::NUM_MILESTONES],
        }
    }

    ///creates a bwt block from the given milestone and bit-vector data.
    pub  (crate)  fn from_data(
        bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
        milestones: [u64; Self::NUM_MILESTONES],
    ) -> Self {
        NucleotideBwtBlock {
            bit_vectors,
            milestones,
        }
    }

    ///Returns the symbol at the given position in the BWT.
    pub  (crate) fn symbol_at(&self, position_block: u64)->Symbol{
        let mut bit_vector_encoding: u64 = 0;

        for bit in 0..self.bit_vectors.len() {
            let bit_value = self.bit_vectors[bit].extract_bit(&position_block);
            bit_vector_encoding |= bit_value << bit;
        }

        Symbol::new_bit_vector(SymbolAlphabet::Nucleotide, bit_vector_encoding as u8)
    }

    ///Sets the symbol at the given position in the BWT block.
    pub (crate) fn set_symbol_at(&mut self, symbol:&Symbol, position_in_block: u64){
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
    #[inline]
    pub (crate) fn set_milestones(&mut self, values: &Vec<u64>) {
        debug_assert!(values.len() >= SymbolAlphabet::Nucleotide.cardinality() as usize);

        unsafe{
            for milestone_idx in 0..SymbolAlphabet::Nucleotide.cardinality() as usize {
                *self.milestones.get_unchecked_mut(milestone_idx) = *values.get_unchecked(milestone_idx);  
            }
        }
    }

    ///gets this block's milestone corresponding to the given symbol.
    #[inline]
    pub (crate) fn milestone(&self, symbol: &Symbol) -> u64 {
        unsafe{
            return *self.milestones.get_unchecked(symbol.index() as usize);  

        }
    }

    /// gets a reference to the milestones array
    pub  (crate)  fn milestones(&self) -> &[u64] {
        &self.milestones
    }
    /// gets a reference to the bit_vectors array

    pub  (crate) fn bit_vectors(&self) -> &[Vec256] {
        &self.bit_vectors
    }

    /// Gets the result of the occurrence function for the local position in this function.
    /// The occurrence function uses the milestone value and the masked occurrenc vector to
    /// determine how many instances of the given character were before this position.
    #[inline]
    pub  (crate) fn global_occurrence(&self, local_query_position: u64, symbol: &Symbol) -> u64 {
        unsafe{
        let milestone_count = self.milestone(&symbol);
        let vec0 = self.bit_vectors.get_unchecked(0).to_simd();
        let vec1 = self.bit_vectors.get_unchecked(1).to_simd();
        let vec2 = self.bit_vectors.get_unchecked(2).to_simd();
        let occurrence_vector = match &symbol.index() {
            1 => simd_and(vec1, vec2),  //A:    0b110
            2 => simd_and(vec0, vec2), //C:    0b101
            3 => simd_and(vec0, vec1), //G:    0b011
            4 => simd_andnot(vec2, simd_andnot(vec0, vec1)), //N:    0b010
            5 => simd_andnot(vec2, simd_andnot(vec1, vec0)), //T:    0b001
            _ => {
                panic!("illegal letter index given in global occurrence function symbol idx given: {}", symbol.index());
            } //assume every other character is an N, since it's illegal to search for a sentinel
        };

        let popcount = masked_popcount(occurrence_vector, local_query_position);

        return milestone_count + popcount as u64;
    }
    }
}

impl AminoBwtBlock {
    pub (crate) const NUM_MILESTONES: usize = 24;
    pub (crate) const NUM_BIT_VECTORS: usize = 5;

    /// create a new bwt block, with data zeroed out
    pub  (crate) fn new() -> Self {
        AminoBwtBlock {
            bit_vectors: [Vec256::new(); Self::NUM_BIT_VECTORS],
            milestones: [0; Self::NUM_MILESTONES],
        }
    }

    /// create a new bwt block from the given data.
    pub  (crate) fn from_data(
        bit_vectors: [Vec256; Self::NUM_BIT_VECTORS],
        milestones: [u64; Self::NUM_MILESTONES],
    ) -> Self {
        AminoBwtBlock {
            milestones,
            bit_vectors,
        }
    }

    ///Gets the symbol at the given position in the BWT block.
    pub  (crate) fn symbol_at(&self, position_block: u64)->Symbol{
        let mut bit_vector_encoding: u64 = 0;

        unsafe{
            for bit in 0..self.bit_vectors.len() {
                let bit_value = self.bit_vectors.get_unchecked(bit).extract_bit(&position_block); 
                bit_vector_encoding |= bit_value << bit;
            }

        }

        Symbol::new_bit_vector(SymbolAlphabet::Amino, bit_vector_encoding as u8)
    }

    ///Sets the symbol at the given position in the BWT block.
    pub  (crate) fn set_symbol_at(&mut self, symbol:&Symbol, position_in_block: u64){
        //create a bitmask, we'll use this to set the bit with an OR operation
        let mut encoded_symbol = symbol.bit_vector();

        //sets the bits in the bit-vectors based on the position and symbol given
        let mut bit_vector_idx = 0;

        unsafe{
            while encoded_symbol != 0 {
                if encoded_symbol & 0x1 == 1 {
                    self.bit_vectors.get_unchecked_mut(bit_vector_idx).set_bit(&position_in_block); 
                }
                encoded_symbol >>= 1;
                bit_vector_idx += 1;
            }
        }
    }


    /// sets the milestones for this block with the values given.
    #[inline]
    pub (crate) fn set_milestones(&mut self, values: &Vec<u64>) {
        debug_assert!(values.len() >= SymbolAlphabet::Amino.cardinality() as usize);
        
        for milestone_idx in 0..SymbolAlphabet::Amino.cardinality() as usize {
            unsafe{
                *self.milestones.get_unchecked_mut(milestone_idx) = *values.get_unchecked(milestone_idx);   
            }
        }
    }

    /// returns the milestone value corresponding to the given symbol
    #[inline]
    pub  (crate) fn milestone(&self, symbol: &Symbol) -> u64 {
        unsafe{
            return *self.milestones.get_unchecked(symbol.index() as usize); 
        }
    }

    /// returns a slice view of the milestones for this block
    pub  (crate) fn milestones(&self) -> &[u64; Self::NUM_MILESTONES] {
        &self.milestones
    }

    /// returns a slice view of the bit_vectors for this block
    pub  (crate) fn bit_vectors(&self) -> &[Vec256; Self::NUM_BIT_VECTORS] {
        &self.bit_vectors
    }

    /// Gets the result of the occurrence function for the local position in this function.
    /// The occurrence function uses the milestone value and the masked occurrenc vector to
    /// determine how many instances of the given character were before this position.
    #[inline]
    pub (crate) fn global_occurrence(&self, local_query_position: SearchPtr, symbol: &Symbol) -> SearchPtr {
        let milestone_count = self.milestone(symbol);

        unsafe{
            let vec0 = self.bit_vectors.get_unchecked(0).to_simd();
            let vec1 = self.bit_vectors.get_unchecked(1).to_simd();
            let vec2 = self.bit_vectors.get_unchecked(2).to_simd();
            let vec3 = self.bit_vectors.get_unchecked(3).to_simd();   
            let vec4 = self.bit_vectors.get_unchecked(4).to_simd();
            let occurrence_vector = match symbol.index() {
                1 => simd_and(vec2,simd_andnot(vec4, vec3)), //A:    0b01100
                2 => simd_andnot(vec3, simd_and(simd_and(vec0, vec1), vec2)), //C:    0b10111
                3 => simd_andnot(vec4, simd_and(vec0, vec1)), //D:    0b00011
                4 => simd_andnot(vec4, simd_and(vec1, vec2)), //E: 0b00110
                5 => simd_andnot(vec0, simd_and(simd_and(vec1, vec2), vec3)), //F:    0b11110
                6 => simd_andnot(vec2, simd_andnot(vec0, vec4)), //G:    0b11010
                7 => simd_andnot(vec2, simd_and(vec0, simd_and(vec1, vec3))), //H: 0b11011
                8 => simd_andnot(vec2, simd_andnot(vec1, vec4)), //I:    0b11001
                9 => simd_andnot(vec1, simd_andnot(vec3, vec4)), //K:    0b10101
                10 => simd_andnot(vec1, simd_andnot(vec0, vec4)), //L:    0b11100
                11 => simd_andnot(vec1, simd_and(vec3, simd_and(vec2, vec0))), //M:    0b11101
                12 => simd_andnot(simd_or(vec0, vec1), simd_andnot(vec2, vec3)), //N:    0b01000
                13 => simd_and(vec3, simd_andnot(vec4, vec0)), //P:    0b01001,
                14 => simd_andnot(simd_or(vec0, vec1), simd_andnot(vec3, vec2)),//Q:    0b00100
                15 => simd_andnot(vec2, simd_andnot(vec3, vec4)), //R:    0b10011
                16 => simd_and(vec1, simd_andnot(vec4, vec3)), //S:    0b01010
                17 => simd_and(vec0, simd_andnot(vec4, vec2)), //T:    0b00101
                18 => simd_andnot(vec3, simd_andnot(vec0, vec4)), //V:    0b10110
                19 => simd_andnot(simd_or(vec1, vec2), simd_andnot(vec3, vec0)), //W:    0b00001
                20 => simd_and(simd_and(vec0, vec1), simd_and(vec2, vec3)), //Ambiguity character X:  0b11111
                21 => simd_andnot(simd_or(vec0, vec2), simd_andnot(vec3, vec1)), //Y:    0b00010
                // 0b00000 is sentinel, but since you can't search for sentinels, it is not included here.
                _ => {
                    panic!("illegal letter index given in global occurrence function");
                }
        };
        
        let popcount = masked_popcount(occurrence_vector, local_query_position);
        
        return milestone_count + popcount as u64;
        }
    }
}

/// enum representing a BWT, either as Nucleotide symbols or Amino symbols
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, MemSize)]
#[serde(untagged)]
pub (crate) enum Bwt {
    Nucleotide(Vec<NucleotideBwtBlock>),
    Amino(Vec<AminoBwtBlock>),
}



impl Bwt {
    pub (crate)  const NUM_SYMBOLS_PER_BLOCK: u64 = 256;

    /// sets a single symbol at the given position in the bwt bit vectors.
    /// this function is meant to be run for every position in the bwt
    /// as a part of BWT data creation
    pub  (crate) fn set_symbol_at(&mut self, bwt_position: &SearchPtr, symbol: &Symbol) {
        //find the block, byte, and bit of the data we're setting
        let bwt_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK;
        unsafe{
            match self{
                Bwt::Nucleotide(vec) =>  vec.get_unchecked_mut(bwt_block_idx as usize).set_symbol_at(symbol, position_in_block),
                Bwt::Amino(vec) => vec.get_unchecked_mut(bwt_block_idx as usize).set_symbol_at(symbol, position_in_block),
            }
        }
    }

    pub(crate) fn num_blocks(bwt_len: u64) -> usize {
        bwt_len.div_ceil(Self::NUM_SYMBOLS_PER_BLOCK as u64) as usize
    }

    /// reconstructs the symbol stored at the given bwt position
    pub  (crate) fn symbol_at(&self, bwt_position: &SearchPtr) -> Symbol {
        //find the block, byte, and bit of the data we're setting
        let position_block_idx = bwt_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let position_in_block = bwt_position % Self::NUM_SYMBOLS_PER_BLOCK;

        unsafe{
            match &self {
                Bwt::Nucleotide(vec) => {
                    let bwt_block = vec.get_unchecked(position_block_idx as usize);
                    bwt_block.symbol_at(position_in_block)
                }

                Bwt::Amino(vec) => {
                    let bwt_block = vec.get_unchecked(position_block_idx as usize);
                    bwt_block.symbol_at(position_in_block)
                }
            }
        }
    }

    ///sets the milestone values based on the given counts array
    pub  (crate) fn set_milestones(&mut self, block_idx: usize, counts: &Vec<u64>) {
        unsafe{
            match self {
                Bwt::Nucleotide(vec) => vec.get_unchecked_mut(block_idx).set_milestones(counts), 
                Bwt::Amino(vec) => vec.get_unchecked_mut(block_idx).set_milestones(counts), 
            }
        }
    }

    /// finds the total occurrence value for the given symbol at the specified global position
    pub (crate) fn global_occurrence(
        &self,
        pointer_global_position: SearchPtr,
        symbol: &Symbol,
    ) -> SearchPtr {
        let block_idx: u64 = pointer_global_position / Self::NUM_SYMBOLS_PER_BLOCK;
        let local_query_position: u64 = pointer_global_position % Self::NUM_SYMBOLS_PER_BLOCK;
        unsafe{

            match self {
                Bwt::Nucleotide(vec) => {
                    vec.get_unchecked(block_idx as usize).global_occurrence(local_query_position, symbol)
                }
                
                Bwt::Amino(vec) => {
                    vec.get_unchecked(block_idx as usize).global_occurrence(local_query_position, symbol)
                }
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
                let expected_value = counts.get(&(position, symbol_idx)).expect("failed to get value from hash table");
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
