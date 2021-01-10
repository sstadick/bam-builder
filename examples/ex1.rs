use bam_builder::{BamBuilder, ReadSortOrder};

fn main() {
    println!("Welcome to example1");
    // Create a builder with all defaults except the read_len is 100
    let mut builder = BamBuilder::new(
        100,
        30,
        "Pair".to_owned(),
        None,
        ReadSortOrder::Unsorted,
        None,
        None,
    );

    // Create a builder for read pair spec
    let records = builder
        .pair_builder()
        .start1(1)
        .start2(1)
        .unmapped1(true)
        .unmapped2(true)
        .bases1("A".repeat(100))
        .bases2("C".repeat(100))
        .build()
        .unwrap();
    println!("{:?}", records);

    // Add the pair to bam builder
    builder.add_pair(records);
    println!("{:?}", builder)
}
