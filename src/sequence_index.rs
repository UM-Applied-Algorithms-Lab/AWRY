use anyhow::Ok;
use libsufr::SequenceFileData;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
struct SequenceMetadata {
    start_position: usize,
    header: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SequenceIndex {
    sequences: Vec<SequenceMetadata>,
}

#[derive(Ord, PartialOrd, Eq, PartialEq)]
pub struct LocalizedSequencePosition {
    sequence_idx: usize,
    local_position: usize,
}

impl LocalizedSequencePosition {
    pub fn new(sequence_idx: usize, local_position: usize) -> Self {
        LocalizedSequencePosition {
            sequence_idx,
            local_position,
        }
    }
    pub fn local_position(&self) -> usize {
        self.local_position
    }
    pub fn sequence_idx(&self) -> usize {
        self.sequence_idx
    }
}

impl SequenceIndex {
    pub fn new() -> Self {
        SequenceIndex {
            sequences: Vec::new(),
        }
    }
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

    pub fn add_sequence(&mut self, header: &str, start_position: usize) {
        self.sequences.push(SequenceMetadata {
            start_position,
            header: header.to_string(),
        })
    }
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

    pub fn serialize<W: std::io::Write>(&self, writer: &mut W) -> Result<(), anyhow::Error> {
        writer.write_all(&(self.sequences.len() as u64).to_le_bytes())?;
        for sequence in self.sequences.iter() {
            writer.write_all(&(sequence.start_position as u64).to_le_bytes())?;
            writer.write_all(&(sequence.header.len() as u64).to_le_bytes())?;
            writer.write_all(sequence.header.as_bytes())?;
        }
        return Ok(());
    }

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
