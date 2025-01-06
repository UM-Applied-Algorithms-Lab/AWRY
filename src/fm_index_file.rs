

use crate::bwt::{AminoBwtBlock, NucleotideBwtBlock};
use crate::compressed_suffix_array::CompressedSuffixArray;
use crate::fm_index::FmIndex;
use crate::kmer_lookup_table::KmerLookupTable;
use crate::sequence_index::SequenceIndex;
use crate::simd_instructions::Vec256;
use crate::{alphabet::SymbolAlphabet, bwt::Bwt};
use std::convert::TryInto;
use std::io::{self, Read};
use std::slice;
use std::{
    io::{Error, Write},
    path::Path,
};

const FM_FILE_LABEL_STRING:&[u8;11] = b"AWRY-Index\n";

impl FmIndex {
    /// Saves them FM-index to disk at the given file path
    /// 
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use awry::alphabet::SymbolAlphabet;
    /// use std::path::Path;
    /// 
    ///
    /// let build_args = FmBuildArgs {
    ///     input_file_src: "test.fasta".into(),
    ///     suffix_array_output_src: None,
    ///     suffix_array_compression_ratio: Some(16),
    ///     lookup_table_kmer_len: None,
    ///     alphabet: SymbolAlphabet::Nucleotide,
    ///     max_query_len: None,
    ///     remove_intermediate_suffix_array_file: true,
    /// };
    /// let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
    /// fm_index.save(&Path::new("test.awry")).expect("unable to save fm index to file");
    /// ``` 
    pub fn save(&self, file_output_src: &Path) -> Result<(), Error> {
        let mut fm_index_file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true) // Overwrite if exists
            .open(file_output_src)?;
        //write the file label string so that the file type can be identified in a text editor
        fm_index_file.write(FM_FILE_LABEL_STRING)?;
        let header = self.generate_file_header();

        for value in header {
            fm_index_file.write_all(&value.to_le_bytes())?;
        }
        //write the BWT to file
        match self.bwt() {
            crate::bwt::Bwt::Nucleotide(vec) => {
                for block in vec.iter() {
                    for bit_vector in block.bit_vectors() {
                        let bit_vector_values = bit_vector.data();
                        for bit_vector_value in bit_vector_values {
                            fm_index_file.write_all(&bit_vector_value.to_le_bytes())?;
                        }
                    }
                    for milestone in block.milestones() {
                        fm_index_file.write_all(&milestone.to_le_bytes())?;
                    }
                }
            }
            crate::bwt::Bwt::Amino(vec) => {
                for block in vec.iter() {
                    for bit_vector in block.bit_vectors() {
                        let bit_vector_values = bit_vector.data();
                        for bit_vector_value in bit_vector_values {
                            fm_index_file.write_all(&bit_vector_value.to_le_bytes())?;
                        }
                    }
                    for milestone in block.milestones() {
                        fm_index_file.write_all(&milestone.to_le_bytes())?;
                    }
                }
            }
        }

        //write the prefix sums
        for prefix_sum in self.prefix_sums() {
            fm_index_file.write_all(&prefix_sum.to_le_bytes())?;
        }

        //write the sampled suffix array
        for suffix_array_value in self.sampled_suffix_array().data() {
            fm_index_file.write_all(&suffix_array_value.to_le_bytes())?;
        }

        //write the kmer lookup table
        debug_assert!(self.kmer_lookup_table().kmer_len() != 0);
        fm_index_file.write(&self.kmer_lookup_table().kmer_len().to_le_bytes())?;
        for range in self.kmer_lookup_table().table(){
            fm_index_file.write_all(&range.start_ptr.to_le_bytes())?;
            fm_index_file.write_all(&range.end_ptr.to_le_bytes())?;
        }

        self.sequence_index().serialize(&mut fm_index_file).expect("unable to serialize sequence index");

        return Ok(());
    }

    ///Loads the fm-index file from the given file path
    /// 
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use awry::alphabet::SymbolAlphabet;
    /// use std::path::Path;
    /// 
    ///
    /// let build_args = FmBuildArgs {
    ///     input_file_src: "test.fasta".into(),
    ///     suffix_array_output_src: None,
    ///     suffix_array_compression_ratio: Some(16),
    ///     lookup_table_kmer_len: None,
    ///     alphabet: SymbolAlphabet::Nucleotide,
    ///     max_query_len: None,
    ///     remove_intermediate_suffix_array_file: true,
    /// };
    /// let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
    /// fm_index.save(&Path::new("test.awry")).expect("unable to save fm index to file");
    /// 
    /// 
    /// let loaded_fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// ``` 
    pub fn load(fm_file_src: &Path) -> Result<FmIndex, Error> {
        let mut fm_index_file = std::fs::OpenOptions::new()
            .write(false)
            .read(true)
            .open(fm_file_src)?;

        //read and check the file label
        let mut file_label_buffer: [u8; FM_FILE_LABEL_STRING.len()] =
            [0; FM_FILE_LABEL_STRING.len()];
        let mut u64_buffer: [u8; 8] = [0; 8];

        //read the label, and check to make sure it matches what we expect. if it doesn't, it's probably not an fm index file.
        fm_index_file.read_exact(&mut file_label_buffer)?;

        //compare the buffer to the expected file label
        let file_label_validated = file_label_buffer.iter().zip(FM_FILE_LABEL_STRING.iter()).all(|(a,b)| a==b);

        if !file_label_validated{
            return Err(std::io::Error::new(io::ErrorKind::InvalidData, 
                "file provided did not start with expected label, and is probably not an fm index file"));
        }

        //read the header. this section may be rewritten and refactored if more than 1 version is supported
        fm_index_file.read_exact(&mut u64_buffer)?;
        let version_number = u64::from_le_bytes(u64_buffer);

        return FmIndex::read_fm_index_by_version_number(& mut fm_index_file, version_number);

    }


    /// generates the file header from the data in the fm-index. This header
    /// may differ on different versions of the index
    fn generate_file_header(&self) -> [u64;4] {
        match self.version_number() {
            _ => {
                let alphabet_idx: u64 = match self.bwt() {
                    Bwt::Nucleotide(_) => 0,
                    Bwt::Amino(_) => 1,
                };

                //Matches version 1
                return [self.version_number(), 
                self.suffix_array_compression_ratio() as u64, 
                self.bwt_len(), 
                alphabet_idx];

            }
        }
    }

    /// Reads the main contents of an fm index file, depending on the found version number.
    fn read_fm_index_by_version_number(fm_index_file: &mut std::fs::File, version_number: u64) -> Result<FmIndex, Error> {
        match version_number{
            _=>{
                let mut u64_buffer: [u8; 8] = [0; 8];
                //currently only version 1 is supported.
                fm_index_file.read_exact(&mut u64_buffer)?;
                let suffix_array_compression_ratio = usize::from_le_bytes(u64_buffer);

                fm_index_file.read_exact(&mut u64_buffer)?;
                let bwt_len =  u64::from_le_bytes(u64_buffer);

                fm_index_file.read_exact(&mut u64_buffer)?;
                let alphabet_idx =  u64::from_le_bytes(u64_buffer);

                let alphabet = match alphabet_idx{
                    0=>{SymbolAlphabet::Nucleotide},
                    1=>{SymbolAlphabet::Amino},
                    _=>panic!("invalid symbol alphabet , did not match any supported alphabet")
                };

                let suffix_array_num_words = CompressedSuffixArray::compressed_word_len(bwt_len as usize, suffix_array_compression_ratio);
                let num_bwt_blocks = Bwt::num_blocks(bwt_len);
                
                let bwt:Bwt = match alphabet{
                    SymbolAlphabet::Nucleotide => {
                        let mut bwt_block_list = vec![NucleotideBwtBlock::new();num_bwt_blocks];
                        for block_idx in 0..bwt_block_list.len(){
                            const BIT_VECTORS_PER_BLOCK:usize = NucleotideBwtBlock::NUM_BIT_VECTORS;
                            const BIT_VECTOR_BYTES_PER_BLOCK:usize = 32*BIT_VECTORS_PER_BLOCK;
                            const NUM_MILESTONES:usize = NucleotideBwtBlock::NUM_MILESTONES;

                            //read the bit vectors for this block
                            let mut bit_vector_buffer:[u8;BIT_VECTOR_BYTES_PER_BLOCK] = [0;BIT_VECTOR_BYTES_PER_BLOCK];
                            fm_index_file.read_exact(&mut bit_vector_buffer)?;
                            let buffer_ptr = bit_vector_buffer.as_ptr();
                            let vector_ptr = buffer_ptr as *const Vec256;
                            let bit_vector_slice = unsafe{
                                slice::from_raw_parts(vector_ptr, BIT_VECTORS_PER_BLOCK)
                            };

                            //read the milestones for this block
                            let mut milestones:[u64;NUM_MILESTONES] = [0; NUM_MILESTONES];
                            for milestone_idx in 0..milestones.len(){
                                fm_index_file.read_exact(&mut u64_buffer)?;
                                milestones[milestone_idx] = u64::from_le_bytes(u64_buffer);
                            }

                            bwt_block_list[block_idx] = NucleotideBwtBlock::from_data( bit_vector_slice.try_into().unwrap(), milestones); 
                        }

                        Bwt::Nucleotide(bwt_block_list) 

                    },
                    SymbolAlphabet::Amino =>{
                        let mut bwt_block_list = vec![AminoBwtBlock::new();num_bwt_blocks];
                        for block_idx in 0..bwt_block_list.len(){
                            const BIT_VECTORS_PER_BLOCK:usize = AminoBwtBlock::NUM_BIT_VECTORS;
                            
                            //read the bit vectors for this block
                            let mut bit_vector_buffer:[u8;32*BIT_VECTORS_PER_BLOCK] = [0;32*BIT_VECTORS_PER_BLOCK];
                            fm_index_file.read_exact(&mut bit_vector_buffer)?;
                            let buffer_ptr = bit_vector_buffer.as_ptr();
                            let vector_ptr = buffer_ptr as *const Vec256;
                            let bit_vector_slice = unsafe{
                                slice::from_raw_parts(vector_ptr, BIT_VECTORS_PER_BLOCK)
                            };
                            
                            //read the milestones for this block
                            let mut milestones:[u64;AminoBwtBlock::NUM_MILESTONES] = [0; AminoBwtBlock::NUM_MILESTONES];
                            for milestone_idx in 0..milestones.len(){
                                fm_index_file.read_exact(&mut u64_buffer)?;
                                milestones[milestone_idx] = u64::from_le_bytes(u64_buffer);
                            }

                            bwt_block_list[block_idx] = AminoBwtBlock::from_data( bit_vector_slice.try_into().unwrap(), milestones);
                        }
                        Bwt::Amino(bwt_block_list)
                    },
                };

                //write the prefix sums
                let mut prefix_sums:Vec<u64> = vec![0; alphabet.cardinality() as usize+1];
                for prefix_sum_idx in 0..prefix_sums.len() {
                    fm_index_file.read_exact(&mut u64_buffer)?;
                    let prefix_sum_value = u64::from_le_bytes(u64_buffer);
                    prefix_sums[prefix_sum_idx] = prefix_sum_value;
                }

                let mut sampled_suffix_array = CompressedSuffixArray::new(bwt_len as usize, suffix_array_compression_ratio);
                        //write the sampled suffix array
                for suffix_array_idx in 0..suffix_array_num_words{
                    fm_index_file.read_exact(&mut u64_buffer)?;
                    let sa_value = u64::from_le_bytes(u64_buffer);
                    sampled_suffix_array.set_word(sa_value, suffix_array_idx);
                }

                let kmer_lookup_table = KmerLookupTable::from_file(fm_index_file, alphabet)?;

                //read in the sequence index
                let sequence_index = SequenceIndex::from_file(fm_index_file).expect("unable to read sequence index from file");
                return Ok(FmIndex::from_elements(bwt, prefix_sums, sampled_suffix_array, kmer_lookup_table, bwt_len, version_number, sequence_index));
            }
        }
    }
}
