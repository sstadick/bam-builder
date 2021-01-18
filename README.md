![Rust](https://github.com/sstadick/bam-builder/workflows/Rust/badge.svg)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/bam-builder)
[![Crates.io](https://img.shields.io/crates/v/bam-builder.svg)](https://crates.io/crates/bam-builder)

# bam-builder

A library for easily building bams for testing.

## Example

```rust
use bam_builder::{bam_order::BamSortOrder, BamBuilder};

fn main() {
    // Create a builder with all defaults except the read_len is 100
    let mut builder = BamBuilder::new(
        100,                    // default read length
        30,                     // default base quality
        "Pair".to_owned(),      // name of sample
        None,                   // optional read group id
        BamSortOrder::Unsorted, // how to sort reads when `.sort` is called
        None,                   // optional sequence dictionary
        None,                   // optional seed used for generating random bases
    );

    // Create a builder for read pair spec
    let records = builder
        .pair_builder()
        .contig(0)               // reads are mapped to tid 0
        .start1(0)               // start pos of read1
        .start2(200)             // start pos of read2
        .unmapped1(false)        // override default of unmapped
        .unmapped2(false)        // override default of unmapped
        .bases1("A".repeat(100)) // override default random bases with "A"s
        .bases2("C".repeat(100)) // override default random bases with "C"s
        .build()                 // inflate the underlying records and set mate info
        .unwrap();

    // Add the pair to bam builder
    builder.add_pair(records);

    // Write records to a file
    let tmp_file = builder.to_tmp().unwrap();
}
```

## Attributions

This library is more a less a direct clone of [fgbio's SamBuilder](https://github.com/fulcrumgenomics/fgbio/blob/master/src/main/scala/com/fulcrumgenomics/testing/SamBuilder.scala), mixed in with some helper methods from [htsjdk](https://github.com/samtools/htsjdk).
