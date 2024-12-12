# AWRY
**A**vx **W**indowed fm-index in **R**ust? **Y**es!

Generates an Fm-Index of a given biological sequence text (Fasta or Fastq file), and implements Locate() and Search() functionalities.

AWRY is a port of a state-of-the-art, fastest in its class FM-index implementation (https://doi.org/10.1186/s13015-021-00204-6). AWRY supports parallelized searching, with parallel_count() and parallel_locate() functions.

## Building an FM-index
to build an fm-index, create an FmBuildArgs struct, and call FmIndex::new()
```rust
let buildArgs =  FmBuildArgs {
    input_file_src: "my_input.fa",              //sets what the input file for the database text will be
    suffix_array_output_src: None,              //will build to a default location
    suffix_array_compression_ratio: None,       // ratio of suffix array compression, 8 by default
    lookup_table_kmer_len: None,                //by default, chooses reasonable table sizes (Dna=13, Amino=5)
    alphabet: SymbolAlphabet::Nucleotide,       //alphabet to build
    max_query_len: None,                        //if set, only sort suffix array up to n positions
    remove_intermediate_suffix_array_file: true,//deletes the suffix array file if true
}

let fm_index = FmIndex::new(&buildArgs);
```

If you only intend to use the count function, you can set the suffix array compression to a high value like 255 to reduce memory usage. 


## Searching for a query
To search for a query, use to count_string and locate_string functions.

```rust
pub fn count_string(&self, query: &String) -> u64 {
    ...
}

/// Finds the locations in the original text of all isntances of the given query.
pub fn locate_string(&self, query: &String) -> Vec<u64> {
    ...
}
```

# Searching for queries in parallel
To find a large number of queries, searching can be parallelized easily with the parallel_count and parallel_locate functions

```rust
pub fn parallel_count(&self, queries: &Vec<String>) -> Vec<u64> {
    ...
}

// Finds the locations for each query in the query list. This function uses rayon's into_par_iter() for parallelism.
pub fn parallel_locate(&self, queries: &Vec<String>) -> Vec<Vec<u64>> {
    ...
}
```