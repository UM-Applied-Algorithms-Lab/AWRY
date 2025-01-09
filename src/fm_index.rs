use std::path::{Path, PathBuf};

use libsufr::sufr_builder::{SufrBuilder, SufrBuilderArgs};
use libsufr::sufr_file::SufrFile;
use libsufr::util::read_sequence_file;
use mem_dbg::MemSize;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    bwt::{AminoBwtBlock, Bwt, NucleotideBwtBlock},
    compressed_suffix_array::CompressedSuffixArray,
    kmer_lookup_table::KmerLookupTable,
    search::{SearchPtr, SearchRange},
    sequence_index::{LocalizedSequencePosition, SequenceIndex},
};

const FM_VERSION_NUMBER: u64 = 1;

// #[derive(Serialize, Deserialize, Debug)]
///Primary FM-index struct
///
/// # Example
/// ```no_run
/// use awry::fm_index::{FmIndex, FmBuildArgs};
/// use awry::alphabet::SymbolAlphabet;
///
/// let build_args = FmBuildArgs {
///     input_file_src: "test.fasta".into(),
///     suffix_array_output_src: None,
///     suffix_array_compression_ratio: None,
///     lookup_table_kmer_len: None,
///     alphabet: SymbolAlphabet::Nucleotide,
///     max_query_len: None,
///     remove_intermediate_suffix_array_file: false,
/// };
/// let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
/// ```
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, MemSize)]
pub struct FmIndex {
    ///BWT-component of the FM-index
    bwt: Bwt,
    ///Prefix Sums, i.e., where in the BWT the given letter starts its range at
    prefix_sums: Vec<u64>,
    /// suffix array used to resolve kmer locations
    sampled_suffix_array: CompressedSuffixArray,
    /// table with precomputed ranges for the first k substring of a query.
    kmer_lookup_table: KmerLookupTable,
    /// length of the bwt, this is also the length of the suffix array (and len of the full text + 1)
    bwt_len: u64,
    /// version number of this struct implementation
    version_number: u64,
    /// sequence index storing the headers, start positions, and lengths of the sequences
    sequence_index: SequenceIndex,
}

const DEFAULT_SUFFIX_ARRAY_FILE_NAME: &str = "sa.sufr";

///Arguments for builing an FM-index
///
/// # Example
/// ```no_run
/// use awry::fm_index::{FmIndex, FmBuildArgs};  
/// use awry::alphabet::SymbolAlphabet;
///
/// let build_args = FmBuildArgs {
///     input_file_src: "test.fasta".into(), ///     
///     suffix_array_output_src: None,
///     suffix_array_compression_ratio: Some(16),
///     lookup_table_kmer_len: None,
///     alphabet: SymbolAlphabet::Nucleotide,
///     max_query_len: None,
///     remove_intermediate_suffix_array_file: true,
/// };
/// let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
#[derive(Serialize, Deserialize, Debug)]
pub struct FmBuildArgs {
    ///file source for the input, either Fasta or Fastq format
    pub input_file_src: PathBuf,
    ///file source to output the intermediate suffix array file.
    pub suffix_array_output_src: Option<PathBuf>,
    ///How much to downsample the suffix array. downsampling increases locate() time but decreases memory usage
    pub suffix_array_compression_ratio: Option<u8>,
    ///Kmer length in the lookup table. Skips the first k search steps at the cost of exponential memory.
    /// If None, uses sensible defaults.
    pub lookup_table_kmer_len: Option<u8>,
    ///alphabet of the input text, and therefore the FM-index.
    pub alphabet: SymbolAlphabet,
    ///Maximum length to allow for searching, or None for unlimited. Setting this slightly speeds up build times,
    /// but may fail if searching for longer queries than specified/
    pub max_query_len: Option<usize>,
    ///If true, will delete the intermediate suffix array file when done.
    /// This file is unnecessary for proper functionality of the FM-index.
    pub remove_intermediate_suffix_array_file: bool,
}

impl FmBuildArgs {
    /// Creates a new FmBuildArgs struct with the given arguments
    pub fn new(
        input_file_src: PathBuf,
        suffix_array_output_src: Option<PathBuf>,
        suffix_array_compression_ratio: Option<u8>,
        lookup_table_kmer_len: Option<u8>,
        alphabet: SymbolAlphabet,
        max_query_len: Option<usize>,
        remove_intermediate_suffix_array_file: bool,
    ) -> Self {
        FmBuildArgs {
            input_file_src: input_file_src,
            suffix_array_output_src,
            suffix_array_compression_ratio,
            lookup_table_kmer_len,
            alphabet,
            max_query_len,
            remove_intermediate_suffix_array_file,
        }
    }
}

impl FmIndex {
    const DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO: u8 = 8;

    /// Construct a new FM-index using the supplied build args.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use awry::alphabet::SymbolAlphabet;
    ///
    /// let build_args = FmBuildArgs {
    ///     input_file_src: "test.fasta".into(),
    ///     suffix_array_output_src: None,
    ///     suffix_array_compression_ratio: None,
    ///     lookup_table_kmer_len: None,
    ///     alphabet: SymbolAlphabet::Nucleotide,
    ///     max_query_len: None,
    ///     remove_intermediate_suffix_array_file: false,
    /// };
    /// let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
    /// ```
    pub fn new(args: &FmBuildArgs) -> Result<Self, anyhow::Error> {
        let suffix_array_src = args
            .suffix_array_output_src
            .clone()
            .unwrap_or(PathBuf::from(DEFAULT_SUFFIX_ARRAY_FILE_NAME));

        let sequence_delimiter = if args.alphabet == SymbolAlphabet::Nucleotide {
            b'N'
        } else {
            b'X'
        };
        let seq_data = read_sequence_file(&args.input_file_src, sequence_delimiter)?;
        let sequence_index = SequenceIndex::from_seq_file_data(&seq_data);

        let sufr_builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: args.max_query_len,
            is_dna: args.alphabet == SymbolAlphabet::Nucleotide,
            allow_ambiguity: true,
            ignore_softmask: true,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            headers: seq_data.headers,
            num_partitions: 1024, //the sufr library suggests 16 partitions
            sequence_delimiter,
            seed_mask: None,
            random_seed: 0,
        };
        let sufr_builder: SufrBuilder<u64> = SufrBuilder::new(sufr_builder_args)?;
        sufr_builder.write(
            &suffix_array_src
                .to_str()
                .expect("unable to parse suffix array src to string"),
        )?;
        //"low_memory" arg on read(), if true, keeps the SA and sequence on-disk, which is not what we usually want.
        let sufr_file: SufrFile<u64> = SufrFile::read(
            &suffix_array_src
                .to_str()
                .expect("unable to parse suffix array src to string"),
            false,
        )?;
        let bwt_len = sufr_file.len_suffixes;
        let sa_compression_ratio = args
            .suffix_array_compression_ratio
            .unwrap_or(Self::DEFAULT_SUFFIX_ARRAY_COMPRESSION_RATIO);
        let mut compressed_suffix_array =
            CompressedSuffixArray::new(bwt_len as usize, sa_compression_ratio as usize);

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
        let mut suffix_array_file_access = sufr_file.suffix_array_file;
        for (suffix_idx, suffix_array_value) in suffix_array_file_access.iter().enumerate() {
            //generate the sampled suffix array
            if suffix_idx % sa_compression_ratio as usize == 0 {
                compressed_suffix_array.set_value(
                    suffix_array_value,
                    suffix_idx / sa_compression_ratio as usize,
                );
            }
            //set the block milestones, if necessary
            if suffix_idx % Bwt::NUM_SYMBOLS_PER_BLOCK as usize == 0 {
                bwt.set_milestones(
                    suffix_idx / Bwt::NUM_SYMBOLS_PER_BLOCK as usize,
                    &letter_counts,
                );
            }
            sufr_file.text[0];
            // get the letter immediately before the suffix array
            let preceeding_letter_ascii = match suffix_array_value {
                0 => '$',
                _ => sufr_file.text[(suffix_array_value - 1) as usize] as char,
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
        let lookup_table_kmer_len = args
            .lookup_table_kmer_len
            .unwrap_or_else(|| KmerLookupTable::default_kmer_len(args.alphabet));

        //generate  the fm index. This must be done before the kmer lookup table can be populated, so for now just use an empty table
        let mut fm_index = FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array: compressed_suffix_array,
            kmer_lookup_table: KmerLookupTable::empty(lookup_table_kmer_len, args.alphabet),
            bwt_len,
            version_number: FM_VERSION_NUMBER,
            sequence_index,
        };

        //now that the fm index has been generated, we can populate the kmer lookup table
        let mut kmer_lookup_table = KmerLookupTable::new(&fm_index, lookup_table_kmer_len);
        kmer_lookup_table.populate_table(&fm_index);
        fm_index.kmer_lookup_table = kmer_lookup_table;

        if args.remove_intermediate_suffix_array_file {
            std::fs::remove_file(Path::new(&suffix_array_src))?;
        }

        return Ok(fm_index);
    }

    /// Creates a new FM-index from its constituent parts
    pub(crate) fn from_elements(
        bwt: Bwt,
        prefix_sums: Vec<u64>,
        sampled_suffix_array: CompressedSuffixArray,
        kmer_lookup_table: KmerLookupTable,
        bwt_len: u64,
        version_number: u64,
        sequence_index: SequenceIndex,
    ) -> FmIndex {
        FmIndex {
            bwt,
            prefix_sums,
            sampled_suffix_array,
            kmer_lookup_table,
            bwt_len,
            version_number,
            sequence_index,
        }
    }

    /// Gets the alphabet the BWT is made from.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use awry::alphabet::SymbolAlphabet;
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let alphabet = fm_index.alphabet();
    /// ```
    pub fn alphabet(&self) -> SymbolAlphabet {
        match self.bwt() {
            Bwt::Nucleotide(_) => SymbolAlphabet::Nucleotide,
            Bwt::Amino(_) => SymbolAlphabet::Amino,
        }
    }

    /// Gets the suffix array compression ratio from the compressed suffix array.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let suffix_array_compression_ratio = fm_index.suffix_array_compression_ratio();
    /// ```
    pub fn suffix_array_compression_ratio(&self) -> usize {
        self.sampled_suffix_array.compression_ratio()
    }

    /// Gets the length of the BWT
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let bwt_len = fm_index.bwt_len();
    /// ```
    pub fn bwt_len(&self) -> u64 {
        self.bwt_len
    }

    /// Gets the Version number for the FM-index.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let version_number = fm_index.version_number();
    /// ```
    pub fn version_number(&self) -> u64 {
        self.version_number
    }

    /// Gets a reference to the Bwt.
    pub(crate) fn bwt(&self) -> &Bwt {
        &self.bwt
    }

    // Gets a reference to the prefix sums.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let prefix_sums = fm_index.prefix_sums();
    /// ```
    pub fn prefix_sums(&self) -> &Vec<u64> {
        &self.prefix_sums
    }

    /// Gets a reference to the compressed suffix array.
    pub(crate) fn sampled_suffix_array(&self) -> &CompressedSuffixArray {
        &self.sampled_suffix_array
    }

    /// Gets a reference to the kmer lookup table.
    pub(crate) fn kmer_lookup_table(&self) -> &KmerLookupTable {
        return &self.kmer_lookup_table;
    }
    /// Gets a reference to the sequence index.
    pub(crate) fn sequence_index(&self) -> &SequenceIndex {
        return &self.sequence_index;
    }

    /// Finds the search range for the given query. This is the heart of the count() and locate() functions.
    pub(crate) fn get_search_range_for_string(&self, query: &str) -> SearchRange {
        if query.len() < self.kmer_lookup_table.kmer_len() as usize {
            let alphabet = self.alphabet();
            let final_query_index =
                Symbol::new_ascii(alphabet, query.chars().last().unwrap()).index();
            let mut search_range =
                SearchRange::new(self, Symbol::new_index(alphabet, final_query_index));
            for query_char in query.chars().rev().skip(1) {
                if !search_range.is_empty() {
                    let query_symbol = Symbol::new_ascii(alphabet, query_char);
                    search_range = self.update_range_with_symbol(search_range, query_symbol);
                }
            }

            return search_range;
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

    // Finds the counts for each query in the query list. This function uses rayon's par_iter() for parallelism.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use rayon::prelude::*;
    /// use std::path::Path;
    /// let queries = vec!["ACGT", "ACGT"];
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let counts = fm_index.parallel_count(queries.par_iter().copied());    
    /// for count in counts{
    ///     println!("count: {}", count);
    /// }
    /// ```
    pub fn parallel_count<'a>(&self, queries: impl ParallelIterator<Item= &'a str>) -> Vec<u64> 
    
    {
        queries
            .into_par_iter()
            .map(|query| self.count_string(&query))
            .collect()
    }

    // Finds the locations for each query in the query list. This function uses rayon's par_iter() for parallelism.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use rayon::prelude::*;
    /// use std::path::Path;
    /// let queries = vec!["ACGT", "ACGT"];
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let query_location_list = fm_index.parallel_locate(queries.par_iter().copied());
    /// for locations_for_query in query_location_list{
    ///     for location in locations_for_query{
    ///         println!("location: {:?}", location);
    ///     }
    /// }
    /// ```
    pub fn parallel_locate<'a>(
        &self,
        queries: impl ParallelIterator<Item= &'a str>,
    ) -> Vec<Vec<LocalizedSequencePosition>>
    {
        queries
            .into_par_iter()
            .map(|query| self.locate_string(query))
            .collect()
    }

    /// Finds the count for the given query.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let count = fm_index.count_string(&String::from("ACGT"));
    /// ```
    pub fn count_string(&self, query: &str) -> u64 {
        self.get_search_range_for_string(query).len()
    }

    /// Finds the locations in the original text of all isntances of the given query.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let locations = fm_index.locate_string(&String::from("ACGT"));
    /// for location in locations{
    ///     println!("location: {:?}", location);
    /// }
    /// ```
    pub fn locate_string(&self, query: &str) -> Vec<LocalizedSequencePosition> {
        let mut string_locations: Vec<LocalizedSequencePosition> = Vec::new();
        let search_range = self.get_search_range_for_string(query);

        // backstep each location until we find a sampled position
        for search_range_idx in search_range.range_iter() {
            let mut num_backsteps_taken: u64 = 0;
            let mut backstep_position = search_range_idx;
            while !self
                .sampled_suffix_array
                .position_is_sampled(backstep_position as usize)
            {
                backstep_position = self.backstep(backstep_position);
                num_backsteps_taken += 1;
            }

            //read the sampled suffix array value, and add the number of backsteps taken to find the real position
            let sequence_position = self.sampled_suffix_array.reconstruct_value(backstep_position as usize).expect("unable to read from the given suffix array position, this is likely an implementation bug or corrupted data.");
            let location = (sequence_position + num_backsteps_taken) % self.bwt_len;

            let localized_position: LocalizedSequencePosition = self
                .sequence_index
                .get_seq_location(location as usize)
                .expect("unable to find sequence location for given position");
            string_locations.push(localized_position);
        }

        string_locations
    }

    /// Perform a single SearchRange updated using a given symbol.
    ///
    /// # Example
    /// ```no_run
    /// use awry::fm_index::{FmIndex, FmBuildArgs};
    /// use awry::alphabet::{Symbol, SymbolAlphabet};
    /// use awry::search::SearchRange;
    /// use std::path::Path;
    ///
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let mut search_range = SearchRange::new(&fm_index, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A'));  
    /// let updated_search_range = fm_index.update_range_with_symbol(search_range, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'G'));
    /// ```
    pub fn update_range_with_symbol(
        &self,
        search_range: SearchRange,
        query_symbol: Symbol,
    ) -> SearchRange {
        unsafe {
            let query_symbol_idx = query_symbol.index() as usize;
            let letter_prefix_sum = self.prefix_sums.get_unchecked(query_symbol_idx);
            let new_start_ptr = letter_prefix_sum
                + self
                    .bwt
                    .global_occurrence(search_range.start_ptr - 1, &query_symbol);
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
    }

    /// Finds the symbol that preceeds the given search pointer, used for finding the most recent sampled SA position.
    pub fn backstep(&self, search_pointer: SearchPtr) -> SearchPtr {
        let symbol = self.bwt.symbol_at(&search_pointer);
        if symbol.is_sentinel() {
            return 0;
        }
        let symbol_prefix_sum = self.prefix_sums[symbol.index() as usize];
        let global_occurrence = self.bwt.global_occurrence(search_pointer, &symbol);
        return symbol_prefix_sum + global_occurrence - 1;
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::HashMap,
        io::Write,
        ops::Range,
        path::{Path, PathBuf},
    };

    use libsufr::sufr_file::SufrFile;
    use rand::{rngs::StdRng, seq::SliceRandom, thread_rng, Rng, SeedableRng};

    use crate::alphabet::SymbolAlphabet;

    use super::{FmBuildArgs, FmIndex};

    fn compare_index_to_reference(
        fm_index: &FmIndex,
        suffix_array_file_src: &PathBuf,
        test_kmer_len: usize,
    ) {
        let mut kmer_map: HashMap<String, Vec<usize>> = HashMap::new();

        let sufr_file = SufrFile::<u64>::read(
            &suffix_array_file_src
                .to_str()
                .expect("unable to parse suffix array src to string"),
            false,
        )
        .expect("Could not read suffix array file");

        for text_position in 0..sufr_file.text_len.saturating_sub(test_kmer_len as u64) as usize {
            let data_slice = &sufr_file.text[text_position..text_position + test_kmer_len];
            let kmer_string = String::from_utf8(data_slice.to_vec())
                .expect("unable to parse kmer ascii values to string");

            kmer_map
                .entry(kmer_string)
                .or_insert_with(Vec::new)
                .push(text_position);
        }

        //check the count and locate functions for correctness against all kmers in our map

        for key in kmer_map.keys() {
            let count = fm_index.count_string(key);
            let kmer_locations = kmer_map.get(key).unwrap();
            assert_eq!(
                count as usize,
                kmer_locations.len(),
                "count returned by fm index did not match count in the testing kmer map"
            );

            let mut locate_result = fm_index.locate_string(key);
            locate_result.sort();
            let mut reference_positions = kmer_locations.clone();
            reference_positions.sort();

            assert_eq!(
                locate_result.len(),
                reference_positions.len(),
                "lengths of fm locate results and testing comparison results did not match"
            );

            for i in 0..locate_result.len() {
                assert_eq!(locate_result[i].local_position() as usize , reference_positions[i], "position did not match in lists between locate results and reference positions");
            }
        }
    }

    #[test]
    fn test_nucleotide_index() -> anyhow::Result<()> {
        const TEST_KMER_LEN: usize = 24;
        const FASTA_SRC: &str = "test_nucleotide.fasta";
        const FASTA_LINE_LEN: u8 = 80;
        const FASTA_SEQ_LEN: usize = 1847;
        let nucleotide_fasta_src: PathBuf = PathBuf::from("test_nucleotide.fasta");
        let suffix_array_file_src: PathBuf = PathBuf::from("test_nucleotide.sufr");
        let fm_index_src: PathBuf = PathBuf::from("test_nucleotide.awry");

        gen_rand_nucleotide_file(&Path::new(FASTA_SRC), FASTA_LINE_LEN, FASTA_SEQ_LEN)
            .expect("unable to generate random nucleotide fasta");

        //create the fm index
        let fm_index = FmIndex::new(&FmBuildArgs {
            input_file_src: nucleotide_fasta_src,
            suffix_array_output_src: Some(suffix_array_file_src.clone()),
            suffix_array_compression_ratio: None,
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: false,
        })
        .expect("unable to build fm index");

        //save the fm index to file
        fm_index
            .save(&Path::new(&fm_index_src))
            .expect("unable to save fm index to file");

        //create map of kmer->Vec<position>
        compare_index_to_reference(&fm_index, &suffix_array_file_src, TEST_KMER_LEN);

        Ok(())
    }

    #[test]
    fn test_amino_index() -> anyhow::Result<()> {
        //I can't find a good API to get protein sequences, so we're just going to randomly generate a protein fasta

        const TEST_KMER_LEN: usize = 8;
        const FASTA_LINE_LEN: u8 = 80;
        const FASTA_LEN: usize = 300;
        let fasta_src = "amino.fasta".to_owned();
        let suffix_array_file_src = "amino.sufr".to_owned();
        let fm_index_src = "amino.awry".to_owned();

        if !Path::new(&fasta_src).exists() {
            gen_rand_amino_file(&Path::new(&fasta_src), FASTA_LINE_LEN, FASTA_LEN)
                .expect("error when generating the amino fasta file");
        }

        //create the fm index
        let fm_index = FmIndex::new(&FmBuildArgs {
            input_file_src: PathBuf::from(fasta_src),
            suffix_array_output_src: Some(PathBuf::from(&suffix_array_file_src)),
            suffix_array_compression_ratio: None,
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Amino,
            max_query_len: None,
            remove_intermediate_suffix_array_file: false,
        })
        .expect("unable to build fm index");

        //save the fm index to file
        fm_index
            .save(&Path::new(&fm_index_src))
            .expect("unable to save fm index to file");

        //create map of kmer->Vec<position>
        compare_index_to_reference(
            &fm_index,
            &PathBuf::from(suffix_array_file_src),
            TEST_KMER_LEN,
        );

        Ok(())
    }

    #[test]
    fn test_fastq_input() -> anyhow::Result<()> {
        const FASTQ_SRC: &str = "test.fastq";
        const SUFFIX_ARRAY_SRC: &str = "test_fastq.sa";
        const FM_INDEX_SRC: &str = "fastq.awry";

        let num_sequences = 30;
        let mut rng = thread_rng();
        let sequence_lengths: Vec<u64> = (0..num_sequences)
            .map(|_| rng.gen_range(5u64..59u64))
            .collect();
        let sequences = generate_random_fastq(FASTQ_SRC, sequence_lengths);

        let fm_index = FmIndex::new(&FmBuildArgs {
            input_file_src: PathBuf::from(FASTQ_SRC),
            suffix_array_output_src: Some(PathBuf::from(SUFFIX_ARRAY_SRC)),
            suffix_array_compression_ratio: None,
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: false,
        })
        .expect("unable to build fm index");

        //save the fm index to file
        fm_index
            .save(&Path::new(&FM_INDEX_SRC.to_owned()))
            .expect("unable to save fm index to file");

        can_find_all_kmers_in_fasq_index(&fm_index, &sequences);

        Ok(())
    }

    fn can_find_all_kmers_in_fasq_index(fm_index: &FmIndex, sequences: &Vec<String>) {
        for sequence in sequences {
            for kmer_start_idx in 0..sequence.len() {
                let query = sequence[kmer_start_idx..sequence.len()].to_owned();
                let kmer_count = fm_index.count_string(&query);
                assert_ne!(
                    kmer_count, 0,
                    "Kmer count returned zero for kmer in the database sequence list"
                );
            }
        }
    }

    fn gen_rand_nucleotide_file(
        fasta_src: &Path,
        line_len: u8,
        sequence_len: usize,
    ) -> anyhow::Result<()> {
        let mut rng = StdRng::seed_from_u64(0);
        let mut fasta_file = std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(fasta_src)
            .expect("unable to open nucleotide fasta file for writing");

        fasta_file.write(">randomly_generated_fasta\n".to_owned().as_bytes())?;

        let mut letters_remaining = sequence_len;
        while letters_remaining >= line_len as usize {
            for _ in 0..line_len {
                fasta_file.write(index_to_nucleotide(rng.gen_range(0..4)).as_bytes())?;
            }
            fasta_file.write("\n".to_owned().as_bytes())?;
            letters_remaining -= line_len as usize;
        }

        //write the last line
        for _ in 0..letters_remaining {
            fasta_file.write(index_to_nucleotide(rng.gen_range(0..4)).as_bytes())?;
        }
        fasta_file.write("\n".to_owned().as_bytes())?;

        Ok(())
    }

    fn gen_rand_amino_file(
        fasta_src: &Path,
        line_len: u8,
        sequence_len: usize,
    ) -> anyhow::Result<()> {
        let mut rng = StdRng::seed_from_u64(999);
        rng.gen_range(0..21);
        let mut fasta_file = std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(fasta_src)
            .expect("unable to open amino fasta file for writing");

        fasta_file.write(">randomly_generated_fasta\n".to_owned().as_bytes())?;

        let mut letters_remaining = sequence_len;
        while letters_remaining >= line_len as usize {
            for _ in 0..line_len {
                fasta_file.write(index_to_amino(rng.gen_range(0..20)).as_bytes())?;
            }
            fasta_file.write("\n".to_owned().as_bytes())?;
            letters_remaining -= line_len as usize;
        }

        //write the last line
        for _ in 0..letters_remaining {
            fasta_file.write(index_to_amino(rng.gen_range(0..20)).as_bytes())?;
        }
        fasta_file.write("\n".to_owned().as_bytes())?;

        Ok(())
    }

    fn index_to_nucleotide(index: u8) -> String {
        match index {
            0 => "A",
            1 => "C",
            2 => "G",
            3 => "T",
            _ => "N",
        }
        .to_owned()
    }

    fn index_to_amino(index: u8) -> String {
        match index {
            0 => "A",
            1 => "C",
            2 => "D",
            3 => "E",
            4 => "F",
            5 => "G",
            6 => "H",
            7 => "I",
            8 => "K",
            9 => "L",
            10 => "M",
            11 => "N",
            12 => "P",
            13 => "Q",
            14 => "R",
            15 => "S",
            16 => "T",
            17 => "V",
            18 => "W",
            19 => "Y",
            _ => "X",
        }
        .to_owned()
    }

    fn generate_random_fastq(output_src: &str, seq_lengths: Vec<u64>) -> Vec<String> {
        let mut sequences: Vec<String> = Vec::new();
        let quality_char_set:Vec<char> = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~".chars().collect();
        let nucleotides: Vec<char> = "AGCT".chars().collect();

        let mut fastq_file = std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(output_src)
            .expect("unable to open fastq file for writing");

        let mut rng = thread_rng();
        for seq_length in seq_lengths {
            let sequence = (0..seq_length)
                .map(|_| nucleotides.choose(&mut rng).unwrap())
                .collect::<String>();
            let quality_string = (0..seq_length)
                .map(|_| quality_char_set.choose(&mut rng).unwrap())
                .collect::<String>();

            let header_string = format!("@dummy-0:{}+\n", seq_length);
            fastq_file
                .write_all(header_string.as_bytes())
                .expect("could not write to fastq file");
            fastq_file
                .write_all(sequence.as_bytes())
                .expect("could not write to fastq file");
            fastq_file
                .write_all("\n".as_bytes())
                .expect("could not write to fastq file");
            fastq_file
                .write_all("+\n".as_bytes())
                .expect("could not write to fastq file");
            fastq_file
                .write_all(quality_string.as_bytes())
                .expect("could not write to fastq file");
            fastq_file
                .write_all("\n".as_bytes())
                .expect("could not write to fastq file");

            sequences.push(sequence);
        }
        fastq_file.flush().expect("unable to flush file");

        return sequences;
    }

    fn gen_multi_sequence_nucleotide_fasta(
        file_src: &Path,
        line_len: u8,
        num_sequences: usize,
        sequence_length_range: Range<usize>,
    ) -> anyhow::Result<Vec<String>> {
        let mut rng = StdRng::seed_from_u64(999);
        let mut fasta_file = std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(file_src)
            .expect("unable to open fastq file for writing");

        let mut sequences: Vec<String> = Vec::new();

        for seq_idx in 0..num_sequences {
            //write a fasta header to file
            fasta_file
                .write(format!(">sequence_{}", seq_idx).as_bytes())
                .expect("could not write to fastq file");

            let sequence_length = thread_rng().gen_range(sequence_length_range.clone());
            let sequence = (0..sequence_length)
                .map(|_| index_to_nucleotide(rng.gen_range(0..4)))
                .collect::<String>();

            for letter_idx in 0..sequence_length {
                //write ascii letter at letter_idx in sequence
                if letter_idx % line_len as usize == 0 {
                    fasta_file.write("\n".to_owned().as_bytes())?;
                }
                fasta_file.write(
                    sequence
                        .chars()
                        .nth(letter_idx)
                        .unwrap()
                        .to_string()
                        .as_bytes(),
                )?;
            }
            fasta_file.write("\n".to_owned().as_bytes())?;

            sequences.push(sequence);
        }
        //todo, return the strings
        Ok(sequences)
    }

    #[test]
    fn multi_sequence_fasta_test() -> anyhow::Result<()> {
        const FASTA_SRC: &str = "test_multi_sequence.fasta";
        const SUFFIX_ARRAY_SRC: &str = "test_multi_sequence.sufr";
        const FM_INDEX_SRC: &str = "test_multi_sequence.awry";
        const FASTA_LINE_LEN: u8 = 80;
        const NUM_SEQUENCES: usize = 24;
        const SEQUENCE_LENGTH_RANGE: Range<usize> = 10..20;

        let sequences = gen_multi_sequence_nucleotide_fasta(
            &Path::new(FASTA_SRC),
            FASTA_LINE_LEN,
            NUM_SEQUENCES,
            SEQUENCE_LENGTH_RANGE,
        )
        .expect("unable to generate random nucleotide fasta");

        //create the fm index
        let fm_index = FmIndex::new(&FmBuildArgs {
            input_file_src: PathBuf::from(FASTA_SRC),
            suffix_array_output_src: Some(PathBuf::from(SUFFIX_ARRAY_SRC)),
            suffix_array_compression_ratio: None,
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: false,
        })
        .expect("unable to build fm index");

        //save the fm index to file
        fm_index
            .save(&Path::new(&FM_INDEX_SRC))
            .expect("unable to save fm index to file");

        //create map of kmer->Vec<position>
        find_kmers_in_reference(&fm_index, sequences);

        Ok(())
    }

    fn find_kmers_in_reference(fm_index: &FmIndex, sequences: Vec<String>) {
        for sequence in sequences {
            for kmer_start_idx in 0..sequence.len() {
                let query = sequence[kmer_start_idx..sequence.len()].to_owned();
                let kmer_count = fm_index.count_string(&query);
                assert_ne!(
                    kmer_count, 0,
                    "Kmer count returned zero for kmer in the database sequence list"
                );
            }
        }
    }
    #[test]
    fn save_load_equality_test() {
        gen_multi_sequence_nucleotide_fasta(Path::new("large.fa"), 80, 40, 10..201)
            .expect("unable to generate random nucleotide fasta");

        let build_args = FmBuildArgs {
            input_file_src: PathBuf::from("large.fa"),
            suffix_array_output_src: None,
            suffix_array_compression_ratio: Some(8),
            lookup_table_kmer_len: None,
            alphabet: SymbolAlphabet::Nucleotide,
            max_query_len: None,
            remove_intermediate_suffix_array_file: true,
        };
        let fm_index = FmIndex::new(&build_args).expect("unable to build fm index");
        fm_index
            .save(Path::new("fm_large.awry"))
            .expect("unable to save fm index");

        let fm_index2 = FmIndex::load(Path::new("fm_large.awry")).expect("unable to load fm index");
        assert_eq!(fm_index.alphabet(), fm_index2.alphabet());
        assert_eq!(fm_index.bwt_len(), fm_index2.bwt_len());
        assert_eq!(fm_index.version_number(), fm_index2.version_number());
        assert_eq!(fm_index.bwt(), fm_index2.bwt());
        assert_eq!(fm_index.prefix_sums(), fm_index2.prefix_sums());
        assert_eq!(
            fm_index.sampled_suffix_array(),
            fm_index2.sampled_suffix_array()
        );
        for i in 0..fm_index.kmer_lookup_table().table().len() {
            assert_eq!(
                fm_index.kmer_lookup_table().table()[i],
                fm_index2.kmer_lookup_table().table()[i]
            );
        }
        assert_eq!(fm_index.kmer_lookup_table(), fm_index2.kmer_lookup_table());
        assert_eq!(
            fm_index.kmer_lookup_table().kmer_len(),
            fm_index2.kmer_lookup_table().kmer_len()
        );
        assert_eq!(fm_index.sequence_index(), fm_index2.sequence_index());
        assert_eq!(fm_index, fm_index2);
    }
}
