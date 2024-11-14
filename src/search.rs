use crate::fm_index::FmIndex;

/// Type representing a position in the BWT
pub type SearchPtr = u64;

/// Represents the range in the BWT that corresponds to a query. A range is valid (corresponds to at least one position) as long as start_ptr <= end_ptr
#[derive(Clone)]
pub struct SearchRange {
    pub start_ptr: SearchPtr,
    pub end_ptr: SearchPtr,
}

impl SearchRange {
    /// Creates a new SearchRange, representing all positions in the BWT
    pub fn new(fm_index: &FmIndex) -> Self {
        SearchRange {
            start_ptr: 0 as SearchPtr,
            end_ptr: fm_index.bwt_len() - 1 as SearchPtr,
        }
    }

    /// Creates a new SearchRange, representing only the first position (the sentinel character)
    pub fn zero() -> Self {
        SearchRange {
            start_ptr: 0,
            end_ptr: 0,
        }
    }

    ///returns true if the search range doesn't represent any elements.
    pub fn is_empty(&self) -> bool {
        return self.start_ptr > self.end_ptr;
    }

    ///gets the number of elements represented by the search range
    pub fn len(&self) -> SearchPtr {
        if self.start_ptr > self.end_ptr {
            0
        } else {
            self.end_ptr - self.start_ptr + 1
        }
    }

    /// Returns an interator over the BWT positions corresponding to this search range
    pub fn range_iter(&self) -> core::ops::Range<SearchPtr> {
        self.start_ptr..(self.end_ptr + 1)
    }
}
