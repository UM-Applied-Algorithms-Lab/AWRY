use anyhow::Ok;
use libsufr::SequenceFileData;
use serde::{Deserialize, Serialize};

/// Struct representing the metadata for a sequence in the sequence index
/// 
/// # Example
/// ```
/// use sufr_bwt::sequence_index::SequenceMetadata;
/// 
/// let sequence_metadata = SequenceMetadata {
///     start_position: 0,
///     header: "header".to_owned(),
/// };
/// ```
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
struct SequenceMetadata {
    start_position: usize,
    header: String,
}
/// Struct representing a sequence index, which maintains a list of all sequence headers and their start positions
/// 
/// # Example
/// ```
/// use sufr_bwt::sequence_index::{SequenceIndex, SequenceMetadata};
/// 
/// let sequence_metadata = SequenceMetadata {
///     start_position: 0,
///     header: "header".to_owned(),
/// };
/// let sequence_index = SequenceIndex::new();
/// sequence_index.add_sequence(&sequence_metadata.header, sequence_metadata.start_position);
/// ```
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
pub struct SequenceIndex {
    sequences: Vec<SequenceMetadata>,
}

/// Struct representing a localized position in a sequence, and the sequence index that contains it
/// 
/// # Example
/// ```
/// use sufr_bwt::sequence_index::LocalizedSequencePosition;
/// 
/// let localized_sequence_position = LocalizedSequencePosition::new(0, 0);
/// ```
pub struct LocalizedSequencePosition {
    sequence_idx: usize,
    local_position: usize,
}

impl LocalizedSequencePosition {
    /// Creates a new LocalizedSequencePosition
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::LocalizedSequencePosition;
    /// 
    /// let localized_sequence_position = LocalizedSequencePosition::new(0, 0);
    /// ``` 
    pub fn new(sequence_idx: usize, local_position: usize) -> Self {
        LocalizedSequencePosition {
            sequence_idx,
            local_position,
        }
    }
    /// Gets the local position component of the LocalizedSequencePosition
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::LocalizedSequencePosition;
    /// 
    /// let localized_sequence_position = LocalizedSequencePosition::new(0, 0);
    /// let local_position = localized_sequence_position.local_position(); 
    /// ```
    pub fn local_position(&self) -> usize {
        self.local_position
    }
    pub fn sequence_idx(&self) -> usize {
        self.sequence_idx
    }
}

impl SequenceIndex {
    /// Creates a new SequenceIndex
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// 
    /// let sequence_index = SequenceIndex::new();
    /// ```
    pub fn new() -> Self {
        SequenceIndex {
            sequences: Vec::new(),
        }
    }
    /// Creates a new SequenceIndex from a SequenceFileData struct (from the libsufr crate)
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// use libsufr::SequenceFileData;
    /// 
    /// let sequence_file_data = SequenceFileData::new();
    /// let sequence_index = SequenceIndex::from_seq_file_data(&sequence_file_data);
    /// ```
    pub fn from_seq_file_data(seq_file_data: &SequenceFileData) -> Self {
        let mut sequence_index: SequenceIndex = SequenceIndex::new();
        for seq_idx in 0..seq_file_data.headers.len() {
            sequence_index.add_sequence(
                &seq_file_data.headers[seq_idx],
                seq_file_data.start_positions[seq_idx],
            );
        }

        return sequence_index;
    }
    /// Adds a sequence to the sequence index
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// 
    /// let sequence_index = SequenceIndex::new();
    /// sequence_index.add_sequence(&String::from("header"), 0);
    /// ```
    pub fn add_sequence(&mut self, header: &str, start_position: usize) {
        self.sequences.push(SequenceMetadata {
            start_position,
            header: header.to_string(),
        })
    }

    /// Gets the LocalizedSequencePosition for the given global position
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// 
    /// let sequence_index = SequenceIndex::new();
    /// sequence_index.add_sequence(&String::from("header"), 0);
    /// let localized_sequence_position = sequence_index.get_seq_location(0).unwrap();
    /// ```
    pub fn get_seq_location(&self, global_position: usize) -> Option<LocalizedSequencePosition> {
        return self.find_sequence_binary_search(global_position, 0, self.sequences.len() - 1);
    }

    fn find_sequence_binary_search(
        &self,
        global_position: usize,
        range_low: usize,
        range_high: usize,
    ) -> Option<LocalizedSequencePosition> {
        if range_low == range_high {
            return Some(LocalizedSequencePosition::new(
                range_low,
                global_position - self.sequences[range_low].start_position,
            ));
        }
        let midpoint = (range_low + range_high) / 2;
        if self.sequences[midpoint].start_position > global_position {
            self.find_sequence_binary_search(global_position, range_low, midpoint)
        } else if self.sequences[midpoint].start_position < global_position {
            self.find_sequence_binary_search(global_position, midpoint, range_high)
        } else {
            Some(LocalizedSequencePosition {
                sequence_idx: midpoint,
                local_position: global_position - self.sequences[midpoint].start_position,
            })
        }
    }

    /// Serializes the sequence index to a file
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// 
    /// let sequence_index = SequenceIndex::new();
    /// sequence_index.add_sequence(&String::from("header"), 0);
    /// let mut sequence_index_file = std::fs::File::create("sequence_index.awry").expect("unable to create sequence index file");
    /// sequence_index.serialize(&mut sequence_index_file).expect("unable to serialize sequence index");
    /// ```
    pub fn serialize<W: std::io::Write>(&self, writer: &mut W) -> Result<(), anyhow::Error> {
        writer.write_all(&(self.sequences.len() as u64).to_le_bytes())?;
        for sequence in self.sequences.iter() {
            writer.write_all(&(sequence.start_position as u64).to_le_bytes())?;
            writer.write_all(&(sequence.header.len() as u64).to_le_bytes())?;
            writer.write_all(sequence.header.as_bytes())?;
        }
        return Ok(());
    }

    /// Deserializes the sequence index from a file
    /// 
    /// # Example
    /// ```
    /// use sufr_bwt::sequence_index::SequenceIndex;
    /// 
    /// let mut sequence_index_file = std::fs::File::open("sequence_index.awry").expect("unable to open sequence index file");
    /// let sequence_index = SequenceIndex::from_file(&mut sequence_index_file).expect("unable to deserialize sequence index");        
    /// ```
    pub fn from_file<R: std::io::Read>(reader: &mut R) -> Result<Self, anyhow::Error> {
        let mut sequence_index: SequenceIndex = SequenceIndex::new();
        //read the number of sequences as a usize
        let mut u64_buffer: [u8; 8] = [0; 8];
        reader.read_exact(&mut u64_buffer)?;
        let num_sequences = u64::from_le_bytes(u64_buffer);

        for _ in 0..num_sequences {
            //read the start position as a usize
            reader.read_exact(&mut u64_buffer)?;
            let start_position = u64::from_le_bytes(u64_buffer) as usize;
            //read the length of the header as a usize
            reader.read_exact(&mut u64_buffer)?;
            let header_len = u64::from_le_bytes(u64_buffer);
            //read the header. since header_len is non a constant, we need to read the header in chunks
            //and construct a string of length header_len
            let mut header_string = String::new();
            let mut remaining_header_len = header_len;
            while remaining_header_len > 0 {
                let mut header_buffer: [u8; 1024] = [0; 1024];
                if remaining_header_len > 1024 {
                    reader.read_exact(&mut header_buffer)?;
                    remaining_header_len -= 1024;
                } else {
                    reader.read_exact(&mut header_buffer[0..remaining_header_len as usize])?;
                    remaining_header_len = 0;
                }
                //append the buffer to header_string
                header_string.push_str(std::str::from_utf8(&header_buffer).unwrap());
            }

            sequence_index.add_sequence(&header_string, start_position);
        }

        return Ok(sequence_index);
    }
}
