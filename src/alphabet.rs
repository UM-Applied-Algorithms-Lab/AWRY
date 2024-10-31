enum BwtNucleotideSymbol {
    Ascii(char),
    Index(u8),
    BitVector(u8),
}

enum BwtAminoSymbol {
    Ascii(u8),
    Index(u8),
    BitVector(u8),
}

pub enum Alphabet {
    Nucleotide,
    Amino,
}

pub fn alphabet_cardinality(alphabet: &Alphabet)-> u8 {
    match alphabet {
        Alphabet::Nucleotide => 6,
        Alphabet::Amino => 22,
    }
}

pub fn ascii_to_index(symbol: &char, alphabet: &Alphabet) -> u8 {
    match alphabet {
        Alphabet::Nucleotide => match symbol {
            '#' | '$' => 0,
            'A' | 'a' => 1,
            'C' | 'c' => 2,
            'G' | 'g' => 3,
            'T' | 't' | 'U' | 'u' => 5,
            _ => 4,
        },
        Alphabet::Amino => match symbol {
            '#' | '$' => 0,
            'A' | 'a' => 1,
            'C' | 'c' => 2,
            'D' | 'd' => 3,
            'E' | 'e' => 4,
            'F' | 'f' => 5,
            'G' | 'g' => 6,
            'H' | 'h' => 7,
            'I' | 'i' => 8,
            'K' | 'k' => 9,
            'L' | 'l' => 10,
            'M' | 'm' => 11,
            'N' | 'n' => 12,
            'P' | 'p' => 13,
            'Q' | 'q' => 14,
            'R' | 'r' => 15,
            'S' | 's' => 16,
            'T' | 't' => 17,
            'V' | 'v' => 18,
            'W' | 'w' => 19,
            'Y' | 'y' => 21,

            _ => 20,
        },
    }
}

pub fn ascii_to_encoded(symbol: &char, alphabet: &Alphabet) -> u8 {
    match alphabet {
        Alphabet::Nucleotide => match symbol {
            '#' | '$' => 0b100,
            'A' | 'a' => 0b110,
            'C' | 'c' => 0b101,
            'G' | 'g' => 0b011,
            'T' | 't' | 'U' | 'u' => 0b010,

            _ => 0b000,
        },
        Alphabet::Amino => match symbol {
            '#' | '$' => 0b00000,
            'A' | 'a' => 0b01100,
            'C' | 'c' => 0b10111,
            'D' | 'd' => 0b00011,
            'E' | 'e' => 0b00110,
            'F' | 'f' => 0b11110,
            'G' | 'g' => 0b11010,
            'H' | 'h' => 0b11011,
            'I' | 'i' => 0b11001,
            'K' | 'k' => 0b10101,
            'L' | 'l' => 0b11100,
            'M' | 'm' => 0b11101,
            'N' | 'n' => 0b01000,
            'P' | 'p' => 0b01001,
            'Q' | 'q' => 0b00100,
            'R' | 'r' => 0b10011,
            'S' | 's' => 0b01010,
            'T' | 't' => 0b00101,
            'V' | 'v' => 0b10110,
            'W' | 'w' => 0b00001,
            'Y' | 'y' => 0b00010,
            _ => 0b11111,
        },
    }
}
