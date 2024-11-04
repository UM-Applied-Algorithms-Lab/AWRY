use crate::{
    alphabet::{alphabet_cardinality, Symbol, SymbolAlphabet},
    bwt::{AminoBwtBlock, NucleotideBwtBlock},
    search::{self, SearchPtr, SearchRange},
};
use bwt;
use bwt::Bwt;

pub struct FmIndex {
    bwt: Bwt,
    prefix_sums: Vec<u64>,
}

impl FmIndex {
    pub fn new(
        input_file_src: &String,
        suffix_array_output_src: &String,
        bwt_alphabet: &SymbolAlphabet,
        max_read_length: &Option<usize>,
        num_threads: &Option<usize>,
    ) -> Result<Self, anyhow::Error> {
        let create_args = sufr::CreateArgs {
            input: input_file_src.clone(),
            num_partitions: 16, //the default in the sufr library
            max_context: max_read_length.clone(),
            threads: num_threads.clone(),
            output: Some(suffix_array_output_src.clone()),
            is_dna: match bwt_alphabet {
                SymbolAlphabet::Nucleotide => true,
                _ => false,
            },
        };
        sufr::create(&create_args);

        let suffix_array: libsufr::SufrFile<u64> =
            libsufr::SufrFile::read(&suffix_array_output_src)?;
        let bwt_length = suffix_array.num_suffixes;

        //find the number of blocks needed (integer ceiling funciton)
        let num_bwt_blocks = bwt_length.div_ceil(bwt::NUM_POSITIONS_PER_BLOCK);
        assert_eq!((100 as usize).div_ceil(3), 34);
        assert_eq!((101 as usize).div_ceil(3), 34);
        assert_eq!((102 as usize).div_ceil(3), 34);
        assert_eq!((103 as usize).div_ceil(3), 35);

        let mut bwt = match bwt_alphabet {
            SymbolAlphabet::Nucleotide => {
                Bwt::Nucleotide(vec![NucleotideBwtBlock::new(); num_bwt_blocks as usize])
            }
            SymbolAlphabet::Amino => {
                Bwt::Amino(vec![AminoBwtBlock::new(); num_bwt_blocks as usize])
            }
        };

        let alphabet_cardinality = alphabet_cardinality(&bwt_alphabet);
        let mut letter_counts = vec![0; alphabet_cardinality as usize];

        for (suffix_idx, position_in_suffix_array) in suffix_array.suffix_array.clone().enumerate()
        {
            //set the block milestones, if necessary
            if position_in_suffix_array % bwt::NUM_POSITIONS_PER_BLOCK == 0 {
                bwt.set_milestones(
                    (position_in_suffix_array / bwt::NUM_POSITIONS_PER_BLOCK) as usize,
                    &letter_counts,
                );
            }

            // get the letter immediately before the suffix array
            let preceeding_letter = match suffix_idx {
                0 => '$',
                _ => suffix_array
                    .string_at((suffix_idx - 1) as usize, Some(1))
                    .chars()
                    .next()
                    .unwrap(),
            };

            let preceeding_letter_idx = ascii_to_index(&preceeding_letter, &bwt_alphabet);
            let encoded_letter = ascii_to_encoded(&preceeding_letter, &bwt_alphabet);
            bwt.set_symbol_at(position_in_suffix_array, encoded_letter);
            letter_counts[preceeding_letter_idx as usize] += 1;
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

        return Ok(FmIndex { bwt, prefix_sums });
    }

    pub fn len(&self) -> usize {
        return self.prefix_sums[self.prefix_sums.len() - 1] as usize;
    }

    pub fn count_string(&self, query: &String) -> usize {
        let mut search_range = SearchRange::new(self);

        for query_char in query.chars().rev() {
            search_range = self.update_range_with_char(search_range, &query_char);
            if search_range.is_empty() {
                return 0;
            }
        }

        return 0;
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

    pub fn backstep(&self, search_pointer: SearchPtr)->SearchPtr{
        let symbol = self.bwt.get_symbol_at(&search_pointer);
        let symbol_prefix_sum = self.prefix_sums[symbol.index() as usize];
        let global_occurrence = self.bwt.global_occurrence(search_pointer, &symbol);
        return symbol_prefix_sum + global_occurrence - 1;
    }
}
