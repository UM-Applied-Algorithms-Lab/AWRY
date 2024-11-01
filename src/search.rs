use crate::fm_index::FmIndex;

pub struct SearchRange {
    pub start_ptr: usize,
    pub end_ptr: usize,
}

impl SearchRange {
    pub fn new(fm_index: &FmIndex) -> Self {
        SearchRange {
            start_ptr: (0),
            end_ptr: (fm_index.len()),
        }
    }

    ///returns true if the search range doesn't represent any elements.
    pub fn is_empty(&self) -> bool {
        return self.start_ptr > self.end_ptr;
    }

    ///gets the number of elements represented by the search range
    pub fn len(&self) -> usize {
        if self.start_ptr > self.end_ptr {
            0
        } else {
            self.end_ptr - self.start_ptr + 1
        }
    }
}
