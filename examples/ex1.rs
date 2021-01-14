use bam_builder::{bam_order::BamSortOrder, BamBuilder};

fn main() {
    println!("Welcome to example1");
    // Create a builder with all defaults except the read_len is 100
    let mut builder = BamBuilder::new(
        100,
        30,
        "Pair".to_owned(),
        None,
        BamSortOrder::Unsorted,
        None,
        None,
    );

    // Create a builder for read pair spec
    let records = builder
        .pair_builder()
        .contig(0)
        .start1(0)
        .start2(200)
        .unmapped1(false)
        .unmapped2(false)
        .bases1("A".repeat(100))
        .bases2("C".repeat(100))
        .build()
        .unwrap();
    println!("{:?}", records);

    // Add the pair to bam builder
    builder.add_pair(records);
    println!("{:?}", builder);
    builder.to_path(std::path::Path::new(&String::from("./test.bam")));
}
