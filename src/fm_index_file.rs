use crate::bwt;
use crate::fm_index::FmIndex;
use crate::{alphabet::SymbolAlphabet, bwt::Bwt};
use std::str;
use std::io::{self, Read};
use std::{
    io::{Error, Write},
    path::Path,
};

const FM_FILE_LABEL_STRING:&[u8;11] = b"AWRY-Index\n";

impl FmIndex {
    pub fn save(&self, file_output_src: &Path) -> Result<(), Error> {
        let mut fm_index_file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true) // Overwrite if exists
            .open(file_output_src)?;

        fm_index_file.write(FM_FILE_LABEL_STRING);
        let header = self.generate_file_header();

        for value in header {
            fm_index_file.write_all(&value.to_le_bytes())?;
        }
        //write the BWT to file
        match self.bwt() {
            crate::bwt::Bwt::Nucleotide(vec) => {
                for block in vec {
                    for milestone in block.milestones() {
                        fm_index_file.write_all(&milestone.to_le_bytes())?;
                    }
                    for bit_vector in block.bit_vectors() {
                        let bit_vector_values = bit_vector.to_u64s();
                        for bit_vector_value in bit_vector_values {
                            fm_index_file.write_all(&bit_vector_value.to_le_bytes())?;
                        }
                    }
                }
            }
            crate::bwt::Bwt::Amino(vec) => {
                for block in vec {
                    for milestone in block.milestones() {
                        fm_index_file.write_all(&milestone.to_le_bytes())?;
                    }
                    for bit_vector in block.bit_vectors() {
                        let bit_vector_values = bit_vector.to_u64s();
                        for bit_vector_value in bit_vector_values {
                            fm_index_file.write_all(&bit_vector_value.to_le_bytes())?;
                        }
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
        return Ok(());
    }

    pub fn load(fm_file_src: &Path) -> Result<FmIndex, Error> {
        let mut fm_index_file = std::fs::OpenOptions::new()
            .write(false)
            .read(true)
            .open(fm_file_src)?;
        let mut file_label_buffer: [u8; FM_FILE_LABEL_STRING.len()] =
            [0; FM_FILE_LABEL_STRING.len()];
        let mut u64_buffer: [u8; 8] = [0; 8];

        //read the label, and check to make sure it matches what we expect. if it doesn't, it's probably not an fm index file.
        fm_index_file.read_exact(&mut file_label_buffer);

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


    ///generates the file header from the data in the fm-index. This header
    /// may differ on different versions of the index
    fn generate_file_header(&self) -> Vec<u64> {
        match self.version_number() {
            _ => {
                let alphabet_idx: u64 = match self.bwt() {
                    Bwt::Nucleotide(_) => 0,
                    Bwt::Amino(_) => 1,
                };

                //Matches version 1
                let mut header: Vec<u64> = vec![0; 8];
                header[0] = self.version_number();
                header[1] = self.suffix_array_compression_ratio();
                header[2] = self.bwt_len();
                header[3] = alphabet_idx;
                //the remaining 32 bytes are left empty for now

                header
            }
        }
    }
    fn read_fm_index_by_version_number(fm_index_file: &mut std::fs::File, version_number: u64) -> Result<FmIndex, Error> {
        match version_number{
            _=>{
                let mut u64_buffer: [u8; 8] = [0; 8];
                //currently only version 1 is supported.
                fm_index_file.read_exact(&mut u64_buffer)?;
                let suffix_array_compression_ratio = u64::from_le_bytes(u64_buffer);

                fm_index_file.read_exact(&mut u64_buffer)?;
                let bwt_len =  u64::from_le_bytes(u64_buffer);

                fm_index_file.read_exact(&mut u64_buffer)?;
                let alphabet_idx =  u64::from_le_bytes(u64_buffer);

                let alphabet = match alphabet_idx{
                    0=>{SymbolAlphabet::Nucleotide},
                    _=>{SymbolAlphabet::Amino}
                };

                let compressed_suffix_array_len = bwt_len / suffix_array_compression_ratio;
                let num_bwt_blocks = bwt_len.div_ceil(bwt::NUM_POSITIONS_PER_BLOCK);

                todo!("read in the bwt, prefix sums, and sampled SA");

                return Ok(FmIndex::from_elements(bwt, prefix_sums, sampled_suffix_array, suffix_array_compression_ratio, bwt_len, version_number));
            }
        }
    }
}