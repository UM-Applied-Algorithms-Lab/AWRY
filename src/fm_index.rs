use std::path::Path;

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    bwt::{AminoBwtBlock, NucleotideBwtBlock},
    compressed_suffix_array::CompressedSuffixArray,
    search::{SearchPtr, SearchRange},
};
use bwt::Bwt;

pub const FM_VERSION_NUMBER: u64 = 1;

pub struct FmIndex {
    bwt: Bwt,
    prefix_sums: Vec<u64>,
    sampled_suffix_array: CompressedSuffixArray,
    suffix_array_compression_ratio: u64,
    bwt_len: u64,
    version_number: u64,
}

const DEFAULT_SUFFIX_ARRAY_FILE_NAME: &str = "sa.sufr";

pub struct FmBuildArgs {
    input_file_src: String,
    suffix_array_output_src: Option<String>,
    suffix_array_compression_ratio: Option<u8>,
    alphabet: SymbolAlphabet,
    max_query_len: Option<usize>,
    threads: Option<usize>,
}

impl FmBuildArgs {
    pub fn new(
        input_file_src: String,
        suffix_array_output_src: Option<String>,
        suffix_array_compression_ratio: Option<u8>,
        alphabet: SymbolAlphabet,
        max_query_len: Option<usize>,
        threads: Option<usize>,
    ) -> Self {
        return FmBuildArgs {
            input_file_src,
            suffix_array_output_src,
            suffix_array_compression_ratio,
            alphabet,
            max_query_len,
            threads,
        };
    }
}

impl FmIndex {
    const DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO: u8 = 8;

    pub fn new(args: &FmBuildArgs) -> Result<Self, anyhow::Error> {
        let suffix_array_src = match args.suffix_array_output_src.clone() {
            Some(src) => src,
            None => DEFAULT_SUFFIX_ARRAY_FILE_NAME.to_owned(),
        };
        let sufr_create_args = sufr::CreateArgs {
            input: args.input_file_src.to_owned(),
            num_partitions: 16, //this is the default, so okay I guess?
            max_query_len: args.max_query_len,
            threads: args.threads,
            output: Some(suffix_array_src.clone()),
            is_dna: args.alphabet == SymbolAlphabet::Nucleotide,
            allow_ambiguity: true,
            ignore_softmask: true,
        };
        sufr::create(&sufr_create_args)?;

        let suffix_array_file: libsufr::SufrFile<u64> =
            libsufr::SufrFile::read(&suffix_array_src)?;
        let bwt_len = suffix_array_file.num_suffixes;
        let sa_compression_ratio = match args.suffix_array_compression_ratio {
            Some(ratio) => ratio,
            None => Self::DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO,
        };
        let mut compressed_suffix_array =
            CompressedSuffixArray::new(bwt_len as usize, sa_compression_ratio as u64);

        //find the number of blocks needed (integer ceiling funciton)
        let num_bwt_blocks = bwt_len.div_ceil(Bwt::NUM_SYMBOLS_PER_BLOCK);

        let mut bwt = match args.alphabet {
            SymbolAlphabet::Nucleotide => {
                Bwt::Nucleotide(vec![NucleotideBwtBlock::new(); num_bwt_blocks as usize])
            }
            SymbolAlphabet::Amino => {
                Bwt::Amino(vec![AminoBwtBlock::new(); num_bwt_blocks as usize])
            }
        };

        let alphabet_cardinality = args.alphabet.cardinality();
        let mut letter_counts = vec![0; alphabet_cardinality as usize];

        let mut suffix_array = suffix_array_file.suffix_array;
        for (suffix_idx, suffix_array_value) in suffix_array.iter().enumerate() {
            //generate the sampled suffix array
            if suffix_idx % sa_compression_ratio as usize == 0 {
                compressed_suffix_array.set_value(
                    suffix_array_value,
                    suffix_idx / sa_compression_ratio as usize,
                );
            }
            //set the block milestones, if necessary
            if suffix_array_value % Bwt::NUM_SYMBOLS_PER_BLOCK == 0 {
                bwt.set_milestones(
                    (suffix_array_value / Bwt::NUM_SYMBOLS_PER_BLOCK) as usize,
                    &letter_counts,
                );
            }
            suffix_array_file.text[0];
            // get the letter immediately before the suffix array
            let preceeding_letter_ascii = match suffix_idx {
                0 => '$',
                _ => suffix_array_file.text[(suffix_array_value - 1) as usize] as char,
            };

            let preceding_symbol =
                Symbol::new_ascii(args.alphabet.clone(), preceeding_letter_ascii);
            let preceding_letter_idx = preceding_symbol.index();
            bwt.set_symbol_at(&(suffix_idx as u64), &preceding_symbol);
            letter_counts[preceding_letter_idx as usize] += 1;
        }

        //generate the prefix sums for the letter counts
        let mut prefix_sums: Vec<u64> = vec![0; alphabet_cardinality as usize + 1];
        let mut accumulated_value = 0;
        for idx in 0..alphabet_cardinality + 1 {
            prefix_sums[idx as usize] = accumulated_value;
            if idx != alphabet_cardinality {
                accumulated_value += letter_counts[idx as usize]
            }
        }

        //remove the temp uncompressed suffix array file, since we no longer need it.
        std::fs::remove_file(Path::new(&suffix_array_src))?;

        return Ok(FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array: compressed_suffix_array,
            suffix_array_compression_ratio: sa_compression_ratio as u64,
            bwt_len,
            version_number: FM_VERSION_NUMBER,
        });
    }

    pub fn from_elements(
        bwt: Bwt,
        prefix_sums: Vec<u64>,
        sampled_suffix_array: CompressedSuffixArray,
        suffix_array_compression_ratio: u64,
        bwt_len: u64,
        version_number: u64,
    ) -> FmIndex {
        FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array,
            suffix_array_compression_ratio,
            bwt_len,
            version_number,
        }
    }

    pub fn suffix_array_compression_ratio(&self) -> u64 {
        self.suffix_array_compression_ratio
    }

    pub fn bwt_len(&self) -> u64 {
        self.bwt_len
    }

    pub fn version_number(&self) -> u64 {
        self.version_number
    }
    pub fn bwt(&self) -> &Bwt {
        &self.bwt
    }
    pub fn prefix_sums(&self) -> &Vec<u64> {
        &self.prefix_sums
    }

    pub fn sampled_suffix_array(&self) -> &CompressedSuffixArray {
        &self.sampled_suffix_array
    }

    pub fn len(&self) -> usize {
        return self.prefix_sums[self.prefix_sums.len() - 1] as usize;
    }

    fn get_search_range_for_string(&self, query: &String) -> SearchRange {
        let mut search_range = SearchRange::new(self);
        for query_char in query.chars().rev() {
            search_range = self.update_range_with_char(search_range, &query_char);
        }

        search_range
    }

    pub fn count_string(&self, query: &String) -> u64 {
        let search_range = self.get_search_range_for_string(query);
        if search_range.start_ptr > search_range.end_ptr {
            0
        } else {
            (search_range.end_ptr - search_range.start_ptr) + 1
        }
    }

    pub fn locate_string(&self, query: &String) -> Vec<u64> {
        let mut string_locations: Vec<u64> = Vec::new();
        let search_range = self.get_search_range_for_string(query);

        for search_range_idx in search_range.range_iter() {
            let mut num_backsteps_taken: u64 = 0;
            let mut backstep_position = search_range_idx;
            while !self
                .sampled_suffix_array
                .position_is_sampled(backstep_position)
            {
                backstep_position = self.backstep(backstep_position);
                num_backsteps_taken += 1;
            }
            let suffix_array_position = backstep_position + num_backsteps_taken;
            let sequence_position = match self.sampled_suffix_array.get_value(suffix_array_position as usize){
                Some(position) => position,
                None => panic!("unable to read from the given suffix array position, this is likely an implementation bug or corrupted data."),
            };
            string_locations.push(sequence_position);
        }

        string_locations
    }

    pub fn update_range_with_char(
        &self,
        search_range: SearchRange,
        query_char: &char,
    ) -> SearchRange {
        let alphabet = match &self.bwt {
            Bwt::Nucleotide(_) => &SymbolAlphabet::Nucleotide,
            Bwt::Amino(_) => &SymbolAlphabet::Amino,
        };
        let query_symbol = Symbol::new_ascii(alphabet.clone(), *query_char);
        let query_symbol_idx = query_symbol.index() as usize;
        let letter_prefix_sum = self.prefix_sums[query_symbol_idx];
        let new_start_ptr = self
            .bwt
            .global_occurrence(search_range.start_ptr, &query_symbol);
        let new_end_ptr = letter_prefix_sum
            + self
                .bwt
                .global_occurrence(search_range.end_ptr, &query_symbol)
            - 1;

        SearchRange {
            start_ptr: new_start_ptr,
            end_ptr: new_end_ptr,
        }
    }

    pub fn backstep(&self, search_pointer: SearchPtr) -> SearchPtr {
        let symbol = self.bwt.get_symbol_at(&search_pointer);
        let symbol_prefix_sum = self.prefix_sums[symbol.index() as usize];
        let global_occurrence = self.bwt.global_occurrence(search_pointer, &symbol);
        return symbol_prefix_sum + global_occurrence - 1;
    }
}
