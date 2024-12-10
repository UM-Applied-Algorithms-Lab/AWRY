use std::fmt::Display;

use serde::{Deserialize, Serialize};

///Describes how a symbol is encoded, either as ASCII, 1-to-N integer index, or strided bit-vector format
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Copy)]
pub(crate) enum SymbolEncoding {
    Ascii(char),
    Index(u8),
    BitVector(u8),
}

///Alphabet from which symbols come from. Any Fm-index, Bwt, etc should come from the same alphabet.
///
///  # Example
/// ```
/// use awry::alphabet::SymbolAlphabet;
///
/// let nucleotide_alphabet = SymbolAlphabet::Nucleotide;
///
/// match nucleotide_alphabet{
///     SymbolAlphabet::Nucleotide => println!("nucleotide"),
///     SymbolAlphabet::Amino => println!("amino"),
/// }
/// ```
///
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Copy)]
pub enum SymbolAlphabet {
    Nucleotide,
    Amino,
}

impl Display for SymbolAlphabet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                SymbolAlphabet::Nucleotide => "nucleotide",
                SymbolAlphabet::Amino => "amino",
            }
        )
    }
}

impl SymbolAlphabet {
    #[allow(dead_code)]
    pub(crate) fn alphabet_id(&self) -> u8 {
        match self {
            SymbolAlphabet::Nucleotide => 0,
            SymbolAlphabet::Amino => 1,
        }
    }
    #[allow(dead_code)]
    pub(crate) fn from_id(id: u8) -> Self {
        match id {
            0 => SymbolAlphabet::Nucleotide,
            1 => SymbolAlphabet::Amino,
            _ => panic!("invalid alphabet id given"),
        }
    }
}

///Implementation of a symbol, from a given alphabet, with a given encoding.
///
/// # Example
/// ```
/// use awry::alphabet::{Symbol, SymbolAlphabet};
///
/// let symbol_from_ascii = Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A');
/// let symbol_from_index = Symbol::new_index(SymbolAlphabet::Amino, 4);//Amino Acid E
/// ```
#[derive(Debug, PartialEq, Eq)]
pub struct Symbol {
    alphabet: SymbolAlphabet,
    encoding: SymbolEncoding,
}

impl Display for Symbol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} code {}", self.alphabet, self.ascii())
    }
}

impl SymbolAlphabet {
    ///Returns the cardinality, or how many different symbols can occur in this alphabet.
    pub(crate) fn cardinality(&self) -> u8 {
        match self {
            SymbolAlphabet::Nucleotide => 6,
            SymbolAlphabet::Amino => 22,
        }
    }
    /// Returns the number of encoding symbols (i.e., non-ambiguity Nucleotide or Amino codes) in the alphabet.
    /// This is used for building and using the kmer lookup table
    pub(crate) fn num_encoding_symbols(&self) -> u8 {
        self.cardinality() - 2
    }
}

impl Symbol {
    ///Creates a new Symbol from a given ascii letter
    ///
    /// # Example
    /// ```
    /// use awry::alphabet::{Symbol, SymbolAlphabet};
    ///
    /// let symbol_from_ascii = Symbol::new_ascii(SymbolAlphabet::Nucleotide, 'A');
    /// ```
    pub fn new_ascii(alphabet: SymbolAlphabet, ascii: char) -> Symbol {
        Symbol {
            alphabet,
            encoding: SymbolEncoding::Ascii(ascii.to_ascii_uppercase()),
        }
    }
    ///Creates a new Symbol from a given index into the alphabet
    ///     
    /// # Example
    /// ```
    /// use awry::alphabet::{Symbol, SymbolAlphabet};
    ///
    /// let symbol_from_index = Symbol::new_index(SymbolAlphabet::Amino, 5);
    /// println!("symbol_from_index: {}", symbol_from_index);
    /// ```
    pub fn new_index(alphabet: SymbolAlphabet, index: u8) -> Symbol {
        debug_assert!(index < alphabet.cardinality());
        Symbol {
            alphabet,
            encoding: SymbolEncoding::Index(index),
        }
    }
    ///Creates a new Symbol

    pub(crate) fn new_bit_vector(alphabet: SymbolAlphabet, bit_vector: u8) -> Symbol {
        Symbol {
            alphabet,
            encoding: SymbolEncoding::BitVector(bit_vector),
        }
    }

    /// gets the ascii representation of the symbol with the given alphabet and encoding
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn ascii(&self) -> char {
        match self.to_ascii().encoding {
            SymbolEncoding::Ascii(c) => c,
            _ => panic!("unable to get ascii encoding, this should not be possible logically."),
        }
    }

    /// gets the index representation of the symbol with the given alphabet and encoding
    #[inline]
    pub(crate) fn index(&self) -> u8 {
        match self.to_index().encoding {
            SymbolEncoding::Index(i) => i,
            _ => panic!("unable to get ascii encoding, this should not be possible logically."),
        }
    }

    /// gets the bit-vector representation of the symbol with the given alphabet and encoding
    #[inline]
    pub(crate) fn bit_vector(&self) -> u8 {
        match self.to_bit_vector().encoding {
            SymbolEncoding::BitVector(b) => b,
            _ => panic!("unable to get ascii encoding, this should not be possible logically."),
        }
    }

    /// converts the symbol into an index-encoded symbol
    pub(crate) fn to_index(&self) -> Symbol {
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
    pub(crate) fn to_bit_vector(&self) -> Symbol {
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
                        0 => 0b100, //sentinel
                        1 => 0b110, //A
                        2 => 0b101, //C
                        3 => 0b011, //G
                        5 => 0b001, //T
                        _ => 0b010, //4 is ambiguity character N,
                    },
                    SymbolEncoding::BitVector(encoding) => encoding,
                }),
            },
        }
    }

    /// converts the symbol into an ascii-encoded symbol
    #[allow(dead_code)]
    pub(crate) fn to_ascii(&self) -> Symbol {
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

    /// returns true if the symbol is a sentinel symbol
    pub(crate) fn is_sentinel(&self) -> bool {
        match self.encoding {
            SymbolEncoding::Ascii(val) => val == '$',
            SymbolEncoding::Index(val) => val == 0,
            SymbolEncoding::BitVector(val) => match self.alphabet {
                SymbolAlphabet::Nucleotide => val == 0b100,
                SymbolAlphabet::Amino => val == 0b00000,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Symbol, SymbolAlphabet};

    #[test]
    fn nucleotide_encoding_transform_test() -> anyhow::Result<()> {
        for ascii_letter in "acgtnACGTN$".chars() {
            let alphabet = SymbolAlphabet::Nucleotide;
            let ascii_symbol = Symbol::new_ascii(alphabet, ascii_letter);
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_bit_vector().to_ascii(),
                "ascii->bit_vector->ascii"
            );
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_bit_vector().to_index().to_ascii(),
                "ascii->bitvector->index->ascii"
            );
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_index().to_ascii(),
                "ascii->index->ascii"
            );
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_index().to_bit_vector().to_ascii(),
                "Ascii->index->bitvector->ascii"
            );
            assert_eq!(ascii_symbol, ascii_symbol.to_ascii(), "ascii->ascii");
        }

        Ok(())
    }

    #[test]
    fn amino_encoding_transform_test() -> anyhow::Result<()> {
        for ascii_letter in "acdefghiklmnpqrstvwxynACDEFGHIKLMNPQRSTVWXY$".chars() {
            let alphabet = SymbolAlphabet::Amino;
            let ascii_symbol = Symbol::new_ascii(alphabet, ascii_letter);
            assert_eq!(ascii_symbol, ascii_symbol.to_bit_vector().to_ascii());
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_bit_vector().to_index().to_ascii()
            );
            assert_eq!(ascii_symbol, ascii_symbol.to_index().to_ascii());
            assert_eq!(
                ascii_symbol,
                ascii_symbol.to_index().to_bit_vector().to_ascii()
            );
            assert_eq!(ascii_symbol, ascii_symbol.to_ascii());
        }

        Ok(())
    }
}
