use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    fm_index::FmIndex,
    search::SearchRange,
};
pub struct KmerLookupTable {
    range_table: Vec<SearchRange>,
    kmer_len: u8,
    alphabet: SymbolAlphabet,
}

impl KmerLookupTable {
    pub fn new(kmer_len: u8, fm_index: FmIndex) -> Self {
        let alphabet = fm_index.alphabet();
        let cardinality = alphabet.cardinality();
        let kmer_table_len = cardinality.pow(kmer_len as u32);

        let mut lookup_table = KmerLookupTable {
            range_table: Vec::new(),
            kmer_len,
            alphabet,
        };
        lookup_table.range_table.reserve(kmer_table_len as usize);

        lookup_table.populate_table(&fm_index);

        return lookup_table;
    }

    fn populate_table(&mut self, fm_index: &FmIndex) {
        let search_range = SearchRange::new(&fm_index);
        self.populate_table_recursive(fm_index, &search_range, 1 as usize, 0 as usize, 1 as usize);
    }

    fn populate_table_recursive(
        &mut self,
        fm_index: &FmIndex,
        search_range: &SearchRange,
        current_kmer_len: usize,
        current_kmer_idx: usize,
        letter_multiplier: usize,
    ) {
        //recursive base case
        if current_kmer_len == self.kmer_len as usize {
            self.range_table[current_kmer_idx] = search_range.clone();
            return;
        }

        let alphabet = fm_index.alphabet();
        let cardinality = alphabet.cardinality();

        for idx in 0..cardinality {
            let new_search_range = fm_index
                .update_range_with_symbol(search_range.clone(), Symbol::new_index(alphabet, idx));
            let new_kmer_idx = current_kmer_idx + (idx as usize * letter_multiplier);
            let new_letter_multiplier = letter_multiplier * cardinality as usize;
            self.populate_table_recursive(
                fm_index,
                &new_search_range,
                current_kmer_len + 1,
                new_kmer_idx,
                new_letter_multiplier,
            );
        }
    }
}
