pub enum SymbolEncoding {
    Ascii(char),
    Index(u8),
    BitVector(u8),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SymbolAlphabet {
    Nucleotide,
    Amino,
}

pub struct Symbol {
    alphabet: SymbolAlphabet,
    encoding: SymbolEncoding,
}
impl SymbolAlphabet {
    pub fn cardinality(&self) -> u8 {
        match self {
            SymbolAlphabet::Nucleotide => 6,
            SymbolAlphabet::Amino => 22,
        }
    }
}

impl Symbol {
    pub fn new_ascii(alphabet: SymbolAlphabet, ascii: char) -> Symbol {
        Symbol {
            alphabet,
            encoding: SymbolEncoding::Ascii(ascii),
        }
    }
    pub fn new_index(alphabet: SymbolAlphabet, index: u8) -> Symbol {
        Symbol {
            alphabet,
            encoding: SymbolEncoding::Index(index),
        }
    }
    pub fn new_bit_vector(alphabet: SymbolAlphabet, bit_vector: u8) -> Symbol {
        Symbol {
            alphabet,
            encoding: SymbolEncoding::BitVector(bit_vector),
        }
    }

    /// gets the ascii representation of the symbol with the given alphabet and encoding
    #[inline]
    pub fn ascii(&self) -> char {
        if let SymbolEncoding::Ascii(encoding) = self.to_index().encoding {
            return encoding;
        } else {
            panic!("unable to get index encoding, this should not be possible logically.");
        }
    }

    /// gets the index representation of the symbol with the given alphabet and encoding
    #[inline]
    pub fn index(&self) -> u8 {
        if let SymbolEncoding::Index(encoding) = self.to_index().encoding {
            return encoding;
        } else {
            panic!("unable to get index encoding, this should not be possible logically.");
        }
    }

    /// gets the bit-vector representation of the symbol with the given alphabet and encoding
    #[inline]
    pub fn bit_vector(&self) -> u8 {
        if let SymbolEncoding::BitVector(encoding) = self.to_bit_vector().encoding {
            return encoding;
        } else {
            panic!("unable to get index encoding, this should not be possible logically.");
        }
    }

    /// converts the symbol into an index-encoded symbol
    pub fn to_index(&self) -> Symbol {
        match self.alphabet {
            SymbolAlphabet::Amino => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::Index(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => match encoding {
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
                        _ => 20, //ambiguity character
                    },
                    SymbolEncoding::Index(encoding) => encoding,
                    SymbolEncoding::BitVector(encoding) => match encoding {
                        0b00000 => 0,  //$
                        0b01100 => 1,  //A
                        0b10111 => 2,  //C
                        0b00011 => 3,  //D
                        0b00110 => 4,  //E
                        0b11110 => 5,  //F
                        0b11010 => 6,  //G
                        0b11011 => 7,  //H
                        0b11001 => 8,  //I
                        0b10101 => 9,  //K
                        0b11100 => 10, //L
                        0b11101 => 11, //M
                        0b01000 => 12, //N
                        0b01001 => 13, //P
                        0b00100 => 14, //Q
                        0b10011 => 15, //R
                        0b01010 => 16, //S
                        0b00101 => 17, //T
                        0b10110 => 18, //V
                        0b00001 => 19, //W
                        0b00010 => 21, //Y
                        _ => 20,       // ambiguity char X
                    },
                }),
            },
            SymbolAlphabet::Nucleotide => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::Index(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => match encoding {
                        '#' | '$' => 0,
                        'A' | 'a' => 1,
                        'C' | 'c' => 2,
                        'G' | 'g' => 3,
                        'T' | 't' | 'U' | 'u' => 5,
                        _ => 4, //ambiguity char N
                    },
                    SymbolEncoding::Index(encoding) => encoding,
                    SymbolEncoding::BitVector(encoding) => match encoding {
                        0b100 => 0, //$
                        0b110 => 1, //A
                        0b101 => 2, //C
                        0b011 => 3, //G
                        0b001 => 5, //T
                        _ => 4,     //N
                    },
                }),
            },
        }
    }

    /// converts the symbol into an bit-vector-encoded symbol
    pub fn to_bit_vector(&self) -> Symbol {
        match self.alphabet {
            SymbolAlphabet::Amino => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::BitVector(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => match encoding {
                        '$' | '#' => 0b00000,
                        'a' | 'A' => 0b01100,
                        'c' | 'C' => 0b10111,
                        'd' | 'D' => 0b00011,
                        'e' | 'E' => 0b00110,
                        'f' | 'F' => 0b11110,
                        'g' | 'G' => 0b11010,
                        'h' | 'H' => 0b11011,
                        'i' | 'I' => 0b11001,
                        'k' | 'K' => 0b10101,
                        'l' | 'L' => 0b11100,
                        'm' | 'M' => 0b11101,
                        'n' | 'N' => 0b01000,
                        'p' | 'P' => 0b01001,
                        'q' | 'Q' => 0b00100,
                        'r' | 'R' => 0b10011,
                        's' | 'S' => 0b01010,
                        't' | 'T' => 0b00101,
                        'v' | 'V' => 0b10110,
                        'w' | 'W' => 0b00001,
                        'y' | 'Y' => 0b00010,
                        _ => 0b11111, //ambiguity char X
                    },
                    SymbolEncoding::Index(encoding) => match encoding {
                        0 => 0b00000,
                        1 => 0b01100,
                        2 => 0b10111,
                        3 => 0b00011,
                        4 => 0b00110,
                        5 => 0b11110,
                        6 => 0b11010,
                        7 => 0b11011,
                        8 => 0b11001,
                        9 => 0b10101,
                        10 => 0b11100,
                        11 => 0b11101,
                        12 => 0b01000,
                        13 => 0b01001,
                        14 => 0b00100,
                        15 => 0b10011,
                        16 => 0b01010,
                        17 => 0b00101,
                        18 => 0b10110,
                        19 => 0b00001,
                        21 => 0b00010,
                        _ => 0b11111, //ambiguity char x
                    },
                    SymbolEncoding::BitVector(encoding) => encoding,
                }),
            },
            SymbolAlphabet::Nucleotide => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::BitVector(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => match encoding {
                        '#' | '$' => 0b100,
                        'A' | 'a' => 0b110,
                        'C' | 'c' => 0b101,
                        'G' | 'g' => 0b011,
                        'T' | 't' | 'U' | 'u' => 0b001,
                        _ => 0b010, //ambiguity char N
                    },
                    SymbolEncoding::Index(encoding) => match encoding {
                        0 => 0b100,
                        1 => 0b110,
                        2 => 0b101,
                        3 => 0b011,
                        5 => 0b010,
                        _ => 0b000, //4 is ambiguity character,
                    },
                    SymbolEncoding::BitVector(encoding) => encoding,
                }),
            },
        }
    }

    /// converts the symbol into an ascii-encoded symbol
    pub fn to_ascii(&self) -> Symbol {
        match self.alphabet {
            SymbolAlphabet::Amino => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::Ascii(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => encoding,
                    SymbolEncoding::Index(encoding) => match encoding {
                        0 => '$',
                        1 => 'A',
                        2 => 'C',
                        3 => 'D',
                        4 => 'E',
                        5 => 'F',
                        6 => 'G',
                        7 => 'H',
                        8 => 'I',
                        9 => 'K',
                        10 => 'L',
                        11 => 'M',
                        12 => 'N',
                        13 => 'P',
                        14 => 'Q',
                        15 => 'R',
                        16 => 'S',
                        17 => 'T',
                        18 => 'V',
                        19 => 'W',
                        21 => 'Y',
                        _ => 'X', //20, is ambiguity
                    },
                    SymbolEncoding::BitVector(encoding) => match encoding {
                        0b00000 => '$',
                        0b01100 => 'A',
                        0b10111 => 'C',
                        0b00011 => 'D',
                        0b00110 => 'E',
                        0b11110 => 'F',
                        0b11010 => 'G',
                        0b11011 => 'H',
                        0b11001 => 'I',
                        0b10101 => 'K',
                        0b11100 => 'L',
                        0b11101 => 'M',
                        0b01000 => 'N',
                        0b01001 => 'P',
                        0b00100 => 'Q',
                        0b10011 => 'R',
                        0b01010 => 'S',
                        0b00101 => 'T',
                        0b10110 => 'V',
                        0b00001 => 'W',
                        0b00010 => 'Y',
                        _ => 'X',
                    },
                }),
            },
            SymbolAlphabet::Nucleotide => Symbol {
                alphabet: self.alphabet,
                encoding: SymbolEncoding::Ascii(match self.encoding {
                    SymbolEncoding::Ascii(encoding) => encoding,
                    SymbolEncoding::Index(encoding) => match encoding {
                        0 => '$',
                        1 => 'A',
                        2 => 'C',
                        3 => 'G',
                        5 => 'T',
                        _ => 'N',
                    },
                    SymbolEncoding::BitVector(encoding) => match encoding {
                        0b100 => '$',
                        0b110 => 'A',
                        0b101 => 'C',
                        0b011 => 'G',
                        0b001 => 'T',
                        _ => 'N',
                    },
                }),
            },
        }
    }
}
