pub mod simd_instructions;
pub mod alphabet;
pub mod bwt;
pub mod fm_index;
pub mod search;
pub mod occurrence;
use std::intrinsics;
use sufr;

//todo: how am I going to square the Suffix Array and the read collection?
//how to get the suffix array data out (use the sufr library its self)
//todo: count function, returns range, method on range to get len
//todo: set up hashmap for cannonical kmers
//todo: allow for adding 0-weighted references (linear walk through fastas)

fn create_sufr_suffix_array()-> Result<()>{

}