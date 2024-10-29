use alphabet::Alphabet;
use std::ops::BitAnd;
use std::ops::BitOr;
use std::ops::BitXor;
use std::simd::Simd;

use crate::alphabet::{alphabet_cardinality, ascii_to_encoded, ascii_to_index};
type SimdVec256 = wide::u64x4;

pub const NUM_NUCLEOTIDE_MILESTONES: usize = 8; //4 nucs, plus N, and $, padded to 8
pub const NUM_AMINO_MILESTONES: usize = 24; //20 AAs, plus X, and $, padded to 24
pub const NUM_NUCLEOTIDE_BIT_VECTORS: usize = 3;
pub const NUM_AMINO_BIT_VECTORS: usize = 5;
pub const NUM_POSITIONS_PER_BLOCK: usize = 256;

#[repr(align(32))]
struct AlignedVectorArray {
    data: [u64; 4],
}
impl AlignedVectorArray {
    fn new() -> Self {
        AlignedVectorArray { data: [0; 4] }
    }
}

#[derive(Clone)]
pub struct NucleotideBwtBlock {
    milestones: [u64; NUM_NUCLEOTIDE_MILESTONES],
    bit_vectors: [SimdVec256; NUM_NUCLEOTIDE_BIT_VECTORS],
}

#[derive(Clone)]
pub struct AminoBwtBlock {
    milestones: [u64; NUM_AMINO_MILESTONES],
    bit_vectors: [SimdVec256; NUM_AMINO_BIT_VECTORS],
}

impl NucleotideBwtBlock {
    pub fn new() -> Self {
        let aligned_empty_vector = AlignedVectorArray::new();
        NucleotideBwtBlock {
            milestones: [0; NUM_NUCLEOTIDE_MILESTONES],
            bit_vectors: [SimdVec256::new(aligned_empty_vector.data); NUM_NUCLEOTIDE_BIT_VECTORS],
        }
    }

    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        for milestone_idx in 0..NUM_NUCLEOTIDE_MILESTONES {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }

    #[inline]
    pub fn get_milestone(&self, letter_idx: u8) -> u64 {
        return self.milestones[letter_idx as usize];
    }

    #[inline]
    pub fn global_occurrence(&self, local_query_position: usize, letter_idx: u8) -> usize {
        let milestone_count = self.get_milestone(letter_idx);
        let vectors = &self.bit_vectors;
        todo!("get a andnot working");
        let occurrence_vector = match letter_idx {
            0 => vectors[1].bitand(vectors[2]), //A:    0b110
            1 => vectors[0].bitand(vectors[2]), //C:    0b101
            2 => vectors[0].bitand(vectors[1]), //G:    0b011
            3 => vectors[0].bitand(vectors[1].bitand(vectors[2])), //T:    0b111
            4 => vectors[0],
            _ => {
                panic!("illegal letter index given in global occurrence function");
            } //assume every other character is an N, since it's illegal to search for a sentinel
        };

        //generate the position bitmask
        let mut position_bitmask_vector = AlignedVectorArray::new();
        let num_fill_bitmask_elements = local_query_position / 64;
        for element_idx in 0..num_fill_bitmask_elements{
            position_bitmask_vector.data[element_idx] = !(0 as u64); //all set ones for a 64-bit int.
        }
        //for the final element, fill with as many ones as needed
        position_bitmask_vector.data[num_fill_bitmask_elements] = !(0 as u64) >> (63 - (local_query_position % 64));
        
        let position_bitmask: SimdVec256 = SimdVec256::new(position_bitmask_vector.data);
        //apply the bitmask and count up the ones in the SIMD vector, giving us the local character count.
        let masked_position_vector = occurrence_vector.bitand(position_bitmask);
        let masked_position_array = masked_position_vector.to_array();

        //get the popcount
        let mut popcount = 0;
        for element_idx in 0..4 {
            popcount += masked_position_array[element_idx].count_ones();
        };

        return popcount as usize + milestone_count as usize;
    }
}


impl AminoBwtBlock {
    pub fn new() -> Self {
        let aligned_empty_vector = AlignedVectorArray::new();
        AminoBwtBlock {
            milestones: [0; NUM_AMINO_MILESTONES],
            bit_vectors: [SimdVec256::new(aligned_empty_vector.data); NUM_AMINO_BIT_VECTORS],
        }
    }
    #[inline]
    pub fn set_milestones(&mut self, values: &Vec<u64>) {
        for milestone_idx in 0..NUM_AMINO_MILESTONES {
            self.milestones[milestone_idx] = values[milestone_idx];
        }
    }
    #[inline]
    pub fn get_milestone(&self, letter_idx: u8) -> u64 {
        return self.milestones[letter_idx as usize];
    }
    #[inline]
    pub fn global_occurrence(&self, local_query_position: usize, letter_idx: u8) -> usize {
        let milestone_count = self.get_milestone(letter_idx);
        let vectors = &self.bit_vectors;
        todo!("get a andnot working");
        todo!("replace occurrence calc ")
        let occurrence_vector = match letter_idx {
            0 => vectors[1].bitand(vectors[2]), //A:    0b110
            1 => vectors[0].bitand(vectors[2]), //C:    0b101
            2 => vectors[0].bitand(vectors[1]), //G:    0b011
            3 => vectors[0].bitand(vectors[1].bitand(vectors[2])), //T:    0b111
            4 => vectors[0],
            _ => {
                panic!("illegal letter index given in global occurrence function");
            } //assume every other character is an N, since it's illegal to search for a sentinel
        };

        //generate the position bitmask
        let mut position_bitmask_vector = AlignedVectorArray::new();
        let num_fill_bitmask_elements = local_query_position / 64;
        for element_idx in 0..num_fill_bitmask_elements{
            position_bitmask_vector.data[element_idx] = !(0 as u64); //all set ones for a 64-bit int.
        }
        //for the final element, fill with as many ones as needed
        position_bitmask_vector.data[num_fill_bitmask_elements] = !(0 as u64) >> (63 - (local_query_position % 64));
        
        let position_bitmask: SimdVec256 = SimdVec256::new(position_bitmask_vector.data);
        //apply the bitmask and count up the ones in the SIMD vector, giving us the local character count.
        let masked_position_vector = occurrence_vector.bitand(position_bitmask);
        let masked_position_array = masked_position_vector.to_array();

        //get the popcount
        let mut popcount = 0;
        for element_idx in 0..4 {
            popcount += masked_position_array[element_idx].count_ones();
        };

        return popcount as usize + milestone_count as usize;
    }
}

pub enum Bwt {
    Nucleotide(Vec<NucleotideBwtBlock>),
    Amino(Vec<AminoBwtBlock>),
}

impl Bwt {
    /// sets a single symbol at the given position in the bwt bit vectors.
    /// this function is meant to be run for every position in the bwt
    /// as a part of BWT data creation
    pub fn set_symbol_at(&self, btw_position: usize, mut encoded_symbol: u8) {
        //find the block, byte, and bit of the data we're setting
        let position_block_idx: usize = (btw_position / NUM_POSITIONS_PER_BLOCK) as usize;
        let position_in_block = btw_position % NUM_POSITIONS_PER_BLOCK;
        //each data element in the SIMD vector is 64 bits, so we're finding the element
        //and bit in the element here.
        let element_position_in_block: usize = position_in_block as usize / 64;
        let bit_position_in_element = position_in_block % 64;

        //create a bitmask, we'll use this to set the bit with an OR operation

        let mut bitmask_array = AlignedVectorArray::new();
        bitmask_array.data[element_position_in_block] = 1 << bit_position_in_element;
        let vector_bitmask = SimdVec256::new(bitmask_array.data);

        let mut bit_vector_idx = 0;
        while encoded_symbol != 0 {
            if encoded_symbol & 0x1 == 0 {
                match self {
                    Bwt::Nucleotide(vec) => {
                        vec[position_block_idx].bit_vectors[bit_vector_idx].bitor(vector_bitmask);
                    }
                    Bwt::Amino(vec) => {
                        vec[position_block_idx].bit_vectors[bit_vector_idx].bitor(vector_bitmask);
                    }
                }
            }
            encoded_symbol >>= 1;
            bit_vector_idx += 1;
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
    pub fn global_occurrence(&self, pointer_global_position: usize, letter_idx: u8) -> usize {
        let block_idx: usize = pointer_global_position / NUM_POSITIONS_PER_BLOCK;
        let local_query_position: usize = pointer_global_position % NUM_POSITIONS_PER_BLOCK;

        match self {
            Bwt::Nucleotide(vec) => {
                vec[block_idx].global_occurrence(local_query_position, letter_idx)
            }
            Bwt::Amino(vec) => vec[block_idx].global_occurrence(local_query_position, letter_idx),
        }
    }
}
