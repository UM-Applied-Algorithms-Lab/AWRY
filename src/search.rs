use crate::fm_index::FmIndex;

/// Type representing a position in the BWT
pub type SearchPtr = u64;

#[derive(Clone)]
pub struct SearchRange {
    pub start_ptr: SearchPtr,
    pub end_ptr: SearchPtr,
}

impl SearchRange {
    pub fn new(fm_index: &FmIndex) -> Self {
        SearchRange {
            start_ptr: (0),
            end_ptr: (fm_index.len() as SearchPtr),
        }
    }

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

    pub fn range_iter(&self) -> core::ops::Range<SearchPtr> {
        self.start_ptr..(self.end_ptr + 1)
    }
}
