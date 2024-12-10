use serde::{Deserialize, Serialize};

use crate::{alphabet::Symbol, fm_index::FmIndex};

/// Type representing a position in the BWT
pub type SearchPtr = u64;

/// Represents the range in the BWT that corresponds to a query. A range is valid (corresponds to at 
/// least one position) as long as start_ptr <= end_ptr
/// 
/// # Example
/// ```
/// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};
/// 
/// //new builds the starting range based on the last symbol in the given query
/// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
/// let search_range = SearchRange::new(&fm_index, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A'));
/// ```
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
pub struct SearchRange {
    pub start_ptr: SearchPtr,
    pub end_ptr: SearchPtr,
}

impl SearchRange {
    /// Creates a new SearchRange, representing all positions in the BWT
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};
    /// 
    /// let fm_index = FmIndex::load(&Path::new("test.awry")).expect("unable to load fm index from file");
    /// let search_range = SearchRange::new(&fm_index, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A'));
    /// ```
    pub fn new(fm_index: &FmIndex, symbol:Symbol) -> Self {
        SearchRange {
            start_ptr: fm_index.prefix_sums()[symbol.index() as usize] as SearchPtr,
            end_ptr: fm_index.prefix_sums()[(symbol.index()+1) as usize]  - 1 as SearchPtr,
        }
    }

    /// Creates a new SearchRange, representing an invalid range (i.e., no elements)
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};    
    /// let search_range = SearchRange::zero();
    /// ```
    pub fn zero() -> Self {
        SearchRange {
            start_ptr: 1,
            end_ptr: 0,
        }
    }

    ///returns true if the search range doesn't represent any elements.
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};
    /// 
    /// let search_range = SearchRange::new(&fm_index, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A'));
    /// let is_empty = search_range.is_empty();
    /// ```
    #[inline]
    pub fn is_empty(&self) -> bool {
        return self.start_ptr > self.end_ptr;
    }

    ///gets the number of elements represented by the search range
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};
    /// 
    /// let search_range = SearchRange::new(&fm_index, Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A'));
    /// let len = search_range.len();
    /// ```
    pub fn len(&self) -> SearchPtr {
        match self.is_empty() {
            true => 0,
            false => self.end_ptr - self.start_ptr + 1,
        }
    }

    /// Returns an interator over the BWT positions corresponding to this search range
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::{FmIndex, FmBuildArgs, SymbolAlphabet};
    /// 
    /// let search_range = fm_index.get_search_range_for_string(&String::from("ACGT"));
    /// for search_range_idx in search_range.range_iter(){
    ///     println!("bwt position in search range: {:?}", search_range_idx);
    /// }
    /// ```
    pub fn range_iter(&self) -> core::ops::Range<SearchPtr> {
        match self.is_empty() {
            true => 0..0,
            false => self.start_ptr..(self.end_ptr + 1),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::search::SearchRange;

    #[test]
    fn search_range_zero_test() {
        assert_eq!(
            SearchRange::zero().len(),
            1,
            "empty search range was not length 1"
        );
    }
    #[test]
    fn search_range_empty_test() {
        assert_eq!(
            SearchRange {
                start_ptr: 1,
                end_ptr: 0
            }
            .len(),
            0,
            "search range with sp > ep did not return length 0"
        );
        assert_eq!(
            SearchRange {
                start_ptr: 999,
                end_ptr: 0
            }
            .len(),
            0,
            "search range with sp > ep did not return length 0"
        );
    }

    #[test]
    fn range_check_search_range() {
        let search_range = SearchRange {
            start_ptr: 999,
            end_ptr: 0,
        };
        let range_iter = search_range.range_iter();
        assert_eq!(
            range_iter.count(),
            0,
            "search range did not match expected for valid range"
        );
    }

    #[test]
    fn range_check_empty_search_range() {
        let search_range = SearchRange {
            start_ptr: 500,
            end_ptr: 499,
        };
        let range_iter = search_range.range_iter();
        assert_eq!(
            range_iter.count(),
            0,
            "count of empty range iterator was not zero"
        );
    }
}
