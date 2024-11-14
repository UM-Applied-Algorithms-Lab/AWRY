use std::path::Path;

use aligned_vec::avec;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    bwt::{AminoBwtBlock, Bwt, NucleotideBwtBlock, SIMD_ALIGNMENT_BYTES},
    compressed_suffix_array::CompressedSuffixArray,
    kmer_lookup_table::KmerLookupTable,
    search::{SearchPtr, SearchRange},
};

pub const FM_VERSION_NUMBER: u64 = 1;

///Primary FM-index struct
pub struct FmIndex {
    ///BWT-component of the FM-index
    bwt: Bwt,
    ///Prefix Sums, i.e., where in the BWT the given letter starts its range at
    prefix_sums: Vec<u64>,
    /// suffix array used to resolve kmer locations
    sampled_suffix_array: CompressedSuffixArray,
    /// table with precomputed ranges for the first k substring of a query.
    kmer_lookup_table: KmerLookupTable,

    bwt_len: u64,
    version_number: u64,
}

const DEFAULT_SUFFIX_ARRAY_FILE_NAME: &str = "sa.sufr";

pub struct FmBuildArgs {
    input_file_src: String,
    suffix_array_output_src: Option<String>,
    suffix_array_compression_ratio: Option<u8>,
    lookup_table_kmer_len: Option<u8>,
    alphabet: SymbolAlphabet,
    max_query_len: Option<usize>,
    threads: Option<usize>,
}

impl FmBuildArgs {
    pub fn new(
        input_file_src: String,
        suffix_array_output_src: Option<String>,
        suffix_array_compression_ratio: Option<u8>,
        lookup_table_kmer_len: Option<u8>,
        alphabet: SymbolAlphabet,
        max_query_len: Option<usize>,
        threads: Option<usize>,
    ) -> Self {
        return FmBuildArgs {
            input_file_src,
            suffix_array_output_src,
            suffix_array_compression_ratio,
            lookup_table_kmer_len,
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

        let suffix_array_file: libsufr::SufrFile<u64> = libsufr::SufrFile::read(&suffix_array_src)?;
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
            SymbolAlphabet::Nucleotide => Bwt::Nucleotide(
                avec![[SIMD_ALIGNMENT_BYTES] | NucleotideBwtBlock::new(); num_bwt_blocks as usize],
            ),
            SymbolAlphabet::Amino => Bwt::Amino(avec![
                [SIMD_ALIGNMENT_BYTES] | AminoBwtBlock::new();
                num_bwt_blocks as usize
            ]),
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

        //create the kmer lookup table
        let lookup_table_kmer_len = match args.lookup_table_kmer_len {
            Some(len) => len,
            None => KmerLookupTable::default_kmer_len(args.alphabet),
        };

        //generate  the fm index. This must be done before the kmer lookup table can be populated, so for now just use an empty table
        let mut fm_index = FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array: compressed_suffix_array,
            kmer_lookup_table: KmerLookupTable::empty(lookup_table_kmer_len, args.alphabet),
            bwt_len,
            version_number: FM_VERSION_NUMBER,
        };

        //now that the fm index has been generated, we can populate the kmer lookup table
        let mut kmer_lookup_table = KmerLookupTable::new(&fm_index, lookup_table_kmer_len);
        kmer_lookup_table.populate_table(&fm_index);
        fm_index.kmer_lookup_table = kmer_lookup_table;

        //remove the temp uncompressed suffix array file, since we no longer need it.
        std::fs::remove_file(Path::new(&suffix_array_src))?;

        return Ok(fm_index);
    }

    /// Creates a new FM-index from its constituent parts
    pub fn from_elements(
        bwt: Bwt,
        prefix_sums: Vec<u64>,
        sampled_suffix_array: CompressedSuffixArray,
        kmer_lookup_table: KmerLookupTable,
        bwt_len: u64,
        version_number: u64,
    ) -> FmIndex {
        FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array,
            kmer_lookup_table,
            bwt_len,
            version_number,
        }
    }

    /// Gets the alphabet the BWT is made from.
    pub fn alphabet(&self) -> SymbolAlphabet {
        match self.bwt() {
            Bwt::Nucleotide(_) => SymbolAlphabet::Nucleotide,
            Bwt::Amino(_) => SymbolAlphabet::Amino,
        }
    }

    /// Gets the suffix array compression ratio from the compressed suffix array.
    pub fn suffix_array_compression_ratio(&self) -> u64 {
        self.sampled_suffix_array.compression_ratio()
    }

    /// Gets the length of the BWT
    pub fn bwt_len(&self) -> u64 {
        self.bwt_len
    }

    /// Gets the Version number for the FM-index.
    pub fn version_number(&self) -> u64 {
        self.version_number
    }

    /// Gets a reference to the Bwt.
    pub fn bwt(&self) -> &Bwt {
        &self.bwt
    }

    // Gets a reference to the prefix sums.
    pub fn prefix_sums(&self) -> &Vec<u64> {
        &self.prefix_sums
    }

    /// Gets a reference to the compressed suffix array.
    pub fn sampled_suffix_array(&self) -> &CompressedSuffixArray {
        &self.sampled_suffix_array
    }

    /// Gets a reference to the kmer lookup table.
    pub fn kmer_lookup_table(&self) -> &KmerLookupTable {
        return &self.kmer_lookup_table;
    }

    /// Finds the search range for the given query. This is the heart of the count() and locate() functions.
    pub fn get_search_range_for_string(&self, query: &String) -> SearchRange {
        if query.len() < self.kmer_lookup_table.kmer_len() as usize {
            let mut search_range = SearchRange::new(self);
            for query_char in query.chars().rev() {
                search_range = self.update_range_with_symbol(
                    search_range,
                    Symbol::new_ascii(self.alphabet(), query_char),
                );
            }

            search_range
        } else {
            let mut search_range = self.kmer_lookup_table.get_range_for_kmer(self, query);
            for query_char in query
                .chars()
                .rev()
                .skip(self.kmer_lookup_table.kmer_len() as usize)
            {
                search_range = self.update_range_with_symbol(
                    search_range,
                    Symbol::new_ascii(self.alphabet(), query_char),
                );
            }

            return search_range;
        }
    }

    // Finds the counts for each query in the query list. This function uses rayon's into_par_iter() for parallelism.
    pub fn parallel_count(&self, queries: &Vec<String>) -> Vec<u64> {
        queries
            .into_par_iter()
            .map(|query| self.count_string(&query))
            .collect()
    }

    // Finds the locations for each query in the query list. This function uses rayon's into_par_iter() for parallelism.
    pub fn parallel_locate(&self, queries: &Vec<String>) -> Vec<Vec<u64>> {
        queries
            .into_par_iter()
            .map(|query| self.locate_string(&query))
            .collect()
    }

    /// Finds the count for the given query.
    pub fn count_string(&self, query: &String) -> u64 {
        let search_range = self.get_search_range_for_string(query);
        if search_range.start_ptr > search_range.end_ptr {
            0
        } else {
            (search_range.end_ptr - search_range.start_ptr) + 1
        }
    }

    /// Finds the locations in the original text of all isntances of the given query.
    pub fn locate_string(&self, query: &String) -> Vec<u64> {
        let mut string_locations: Vec<u64> = Vec::new();
        let search_range = self.get_search_range_for_string(query);

        // backstep each location until we find a sampled position
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

            //read the sampled suffix array value, and add the number of backsteps taken to find the real position
            let sequence_position = match self.sampled_suffix_array.get_value(backstep_position as usize){
                Some(position) => position,
                None => panic!("unable to read from the given suffix array position, this is likely an implementation bug or corrupted data."),
            };
            string_locations.push(sequence_position + num_backsteps_taken);
        }

        string_locations
    }

    /// Perform a single SearchRange updated using a given symbol.
    pub fn update_range_with_symbol(
        &self,
        search_range: SearchRange,
        query_symbol: Symbol,
    ) -> SearchRange {
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

    /// Finds the symbol that preceeds the given search pointer, used for finding the most recent sampled SA position.
    pub fn backstep(&self, search_pointer: SearchPtr) -> SearchPtr {
        let symbol = self.bwt.get_symbol_at(&search_pointer);
        let symbol_prefix_sum = self.prefix_sums[symbol.index() as usize];
        let global_occurrence = self.bwt.global_occurrence(search_pointer, &symbol);
        return symbol_prefix_sum + global_occurrence - 1;
    }
}
