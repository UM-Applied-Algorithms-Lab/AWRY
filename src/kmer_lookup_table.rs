use std::{
    fs::File,
    io::{Error, Read},
};

use crate::{
    alphabet::{Symbol, SymbolAlphabet},
    fm_index::FmIndex,
    search::SearchRange,
};

///Table storing precomputed SearchRanges for all kmers of a given length.
#[derive(Clone)]
pub struct KmerLookupTable {
    range_table: Vec<SearchRange>,
    kmer_len: u8,
}

impl KmerLookupTable {
    const DEFAULT_KMER_LEN_NUCLEOTIDE: u8 = 10;
    const DEFAULT_KMER_LEN_AMINO: u8 = 4;

    /// Creates and populates a new table for all kmers of the given length
    pub fn new(fm_index: &FmIndex, kmer_len: u8) -> Self {
        let alphabet = fm_index.alphabet();
        let kmer_table_len = Self::get_num_table_entries(kmer_len, alphabet);

        let mut lookup_table = KmerLookupTable {
            range_table: vec![SearchRange::zero(); kmer_table_len],
            kmer_len,
        };

        lookup_table.populate_table(fm_index);

        return lookup_table;
    }

    /// Creates an empty table with capacity to hold all kmers of the given length. The populate_table() func can later be called to fill in the table
    pub fn empty(kmer_len: u8, alphabet: SymbolAlphabet) -> Self {
        let mut table = KmerLookupTable {
            range_table: Vec::new(),
            kmer_len: 0,
        };
        let table_size = Self::get_num_table_entries(kmer_len, alphabet);
        table.range_table.reserve(table_size);

        return table;
    }

    /// Reads the kmer table from a file pointer. The file must be pointing to the beginning of the table.
    pub fn from_file(file: &mut File, alphabet: SymbolAlphabet) -> Result<Self, Error> {
        let mut u8_buffer: [u8; 1] = [0; 1];
        file.read_exact(&mut u8_buffer)?;
        let kmer_len = u8::from_le_bytes(u8_buffer);

        let kmer_table_num_entries = Self::get_num_table_entries(kmer_len, alphabet);
        let mut table = KmerLookupTable {
            range_table: Vec::new(),
            kmer_len,
        };

        table.range_table.reserve(kmer_table_num_entries);
        let mut u64_buffer: [u8; 8] = [0; 8];
        for table_idx in 0..kmer_table_num_entries {
            file.read_exact(&mut u64_buffer)?;
            table.range_table[table_idx].start_ptr = u64::from_le_bytes(u64_buffer);
            table.range_table[table_idx].end_ptr = u64::from_le_bytes(u64_buffer);
        }

        return Ok(table);
    }

    /// Gets a reference to the table data
    pub fn table(&self) -> &Vec<SearchRange> {
        return &self.range_table;
    }

    /// Gets the length of the kmers referenced in this table
    pub fn kmer_len(&self) -> u8 {
        self.kmer_len
    }

    /// Gets the SearchRange corresponding to the ginve kmer
    pub fn get_range_for_kmer(&self, fm_index: &FmIndex, kmer: &str) -> SearchRange {
        debug_assert!(kmer.len() >= self.kmer_len as usize);

        let alphabet = fm_index.alphabet();
        let mut search_range = SearchRange::new(
            fm_index,
            Symbol::new_ascii(alphabet, kmer.chars().last().unwrap()),
        );
        for char in kmer
            .chars()
            .rev()
            .skip(1)
            .take((self.kmer_len - 1) as usize)
        {
            search_range =
                fm_index.update_range_with_symbol(search_range, Symbol::new_ascii(alphabet, char));
        }

        return search_range;
    }

    /// Computes the number of entries the table needs for the given kmer length and alphabet.
    fn get_num_table_entries(kmer_len: u8, alphabet: SymbolAlphabet) -> usize {
        let cardinality = alphabet.num_encoding_symbols() as usize;
        let num_entries = cardinality.pow(kmer_len as u32);

        return num_entries;
    }

    /// Populates the table with search ranges. It is recommended to reserve the capacity of the table before calling this method.
    pub fn populate_table(&mut self, fm_index: &FmIndex) {
        for symbol_index in 1..fm_index.alphabet().num_encoding_symbols() {
            let symbol = Symbol::new_index(fm_index.alphabet(), symbol_index);
            let search_range = SearchRange::new(fm_index, symbol);
            self.populate_table_recursive(
                fm_index,
                &search_range,
                1,
                symbol_index as usize,
                fm_index.alphabet().num_encoding_symbols() as usize,
            );
        }
        // self.populate_table_recursive(fm_index, &search_range, 1 as usize, 0 as usize, 1 as usize);
    }

    /// Recursive component of populating tables.
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
        let encoding_symbol_idx_range = 1..alphabet.num_encoding_symbols();

        for idx in encoding_symbol_idx_range {
            let new_search_range = fm_index
                .update_range_with_symbol(search_range.clone(), Symbol::new_index(alphabet, idx));
            let new_kmer_idx = current_kmer_idx + (idx as usize * letter_multiplier);
            let new_letter_multiplier =
                letter_multiplier * alphabet.num_encoding_symbols() as usize;
            self.populate_table_recursive(
                fm_index,
                &new_search_range,
                current_kmer_len + 1,
                new_kmer_idx,
                new_letter_multiplier,
            );
        }
    }

    /// Returns the default kmer lengths for each alphabet type.
    pub fn default_kmer_len(alphabet: SymbolAlphabet) -> u8 {
        match alphabet {
            SymbolAlphabet::Nucleotide => Self::DEFAULT_KMER_LEN_NUCLEOTIDE,
            SymbolAlphabet::Amino => Self::DEFAULT_KMER_LEN_AMINO,
        }
    }
}
