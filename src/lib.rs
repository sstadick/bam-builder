#[macro_use]
extern crate derive_builder;
use rust_htslib::bam::{header::{Header, HeaderRecord}, record::Record, cigar::{Cigar, CigarString}, HeaderView};
use std::convert::{TryFrom, TryInto};
use rand::{prelude::*};
use std::collections::HashMap;
use std::rc::Rc;

const DEFAULT_SEED: usize = 42;
const DEFAULT_BASE_QUALITY: u8 = 30;
const DEFAULT_READ_LENGTH: usize = 150;
const DEFAULT_BASES: [char; 4] = ['A', 'C', 'T', 'G'];


/// Wrapper type for a ReadGroupID
pub struct ReadGroupId(String);

impl Default for ReadGroupId {
    fn default() -> Self {
        Self(String::from("A"))
    }
}

/// Wrapper type for a Sample Name
pub struct SampleName(String);

impl Default for SampleName {
    fn default() -> Self {
        Self(String::from("Sample1"))
    }
}

/// Hold relevant info about the sequences, translates into a header
pub struct SequenceMetadataStore {
    sequences: Vec<SequenceMetadata>
}

/// Holds relevant information about reference sequences
// TODO: expand to full sequence dictionary like object
pub struct SequenceMetadata {
    name: String,
    length: usize,
}


impl From<SequenceMetadata> for HeaderRecord<'_> {

    fn from(meta: SequenceMetadata) -> Self {
        let mut rec = HeaderRecord::new("SQ".as_bytes());
        rec.push_tag("SN".as_bytes(), &meta.name);
        rec.push_tag("LN".as_bytes(), &meta.length);
        rec
    }
}

impl Default for SequenceMetadataStore {
    fn default() -> Self {
        let sequences = (1..22).map(|chr| chr.to_string()).chain(vec!["X", "Y", "M"].into_iter().map(|chr| chr.to_string())).map(|chr| {
            SequenceMetadata {
                name: format!("chr{}", chr),
                length: 2_000_000
            }
        }).collect();
        SequenceMetadataStore {
            sequences
        }
    }
}

impl From<SequenceMetadataStore> for Header {
    fn from(meta: SequenceMetadataStore) -> Header {
        let mut header = Header::new();
        for seq in meta.sequences.into_iter() {
            header.push_record(&seq.into());
        }
        header
    }
}


/// Sort order options for reads
pub enum ReadSortOrder {
    NameSorted,
    CoordSorted,
    Unsorted,
}

impl Default for ReadSortOrder {
    fn default() -> Self {
        Self::Unsorted
    }
}

#[derive(Clone)]
pub enum Strand {
    Plus,
    Minus
}

/// Builder for SAM records
// See: https://github.com/fulcrumgenomics/fgbio/blob/57988d615fa6188e0028faedfdf65ed13cf645e7/src/main/scala/com/fulcrumgenomics/testing/SamBuilder.scala
pub struct SamBuilder {
    /// Length of reads to generate
    read_length: usize,
    /// Default base quality to use across all qualities
    base_quality: u8,
    /// Name of the sample
    sample_name: SampleName,
    /// The read group to include in tags
    read_group_id: ReadGroupId,
    /// The sort order of the generated reads
    sort: ReadSortOrder,
    /// The BAM header
    header: Header,
    /// Random number generator
    rng: StdRng,
    /// The collection of records being accumulated
    records: Vec<Record>,
    /// The default string o qualities to use
    default_quals: String,
    /// Counter for helping to increment the name
    counter: usize
}

impl SamBuilder {
    pub fn new(
        read_length: usize,
        base_quality: u8,
        sample_name: String,
        read_group_id: Option<String>,
        sort: ReadSortOrder,
        pseudo_sd: Option<SequenceMetadataStore>,
        seed: Option<usize>,
    ) -> SamBuilder {

        let pseudo_sd = pseudo_sd.unwrap_or(SequenceMetadataStore::default());
        let seed = seed.unwrap_or(DEFAULT_SEED);
        let rng = StdRng::seed_from_u64(seed as u64);
        let read_group_id = match read_group_id {
            Some(rgi) => ReadGroupId(rgi),
            None => ReadGroupId::default()
        };
        let default_quals = std::iter::repeat(char::from(33 + base_quality)).take(read_length).collect::<String>();

        let header: Header  = pseudo_sd.into();

        // TODO: Add read group to header


        SamBuilder {
            read_length,
            base_quality,
            sample_name: SampleName(sample_name),
            read_group_id,
            sort,
            header,
            rng,
            records: vec![],
            default_quals,
            counter: 0,
        }
    }

    /// Generate a sequential name, 16 chars wide padded with 0's
    fn next_name(&mut self) -> String {
        let name = format!("{:0>16}", self.counter);
        self.counter += 1;
        name
    }

    /// Generate a random sequence of bases of length readLength
    fn random_bases(&mut self) -> String {
        (0..self.read_length).map(|_| DEFAULT_BASES[self.rng.gen_range(0, DEFAULT_BASES.len())]).collect::<String>()
    }

    /// Create a ReadPairSpecBuilder with defaults filled in based on the SamBuilder
    pub fn pair_builder(&mut self) -> ReadPairSpecBuilder {
        ReadPairSpecBuilder::default()
            .name(self.next_name())
            .bases1(self.random_bases())
            .bases2(self.random_bases())
            .quals1(self.default_quals.clone())
            .quals2(self.default_quals.clone())
            .contig(0)
            .start1(0) // corresponds to unset TODO: Check this
            .start2(0)
            .unmapped1(true)
            .unmapped2(true)
            .cigar1(format!("{}M", &self.read_length))
            .cigar2(format!("{}M", &self.read_length))
            .mapq1(60)
            .mapq2(60)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .to_owned()
    }

    // Look into using https://stackoverflow.com/questions/19650265/is-there-a-faster-shorter-way-to-initialize-variables-in-a-rust-struct or
    // https://docs.rs/derive_builder/0.5.1/derive_builder/ to collect all the fields to build a bam, and may for the builder itself
    /// Add a read pair.
    ///
    /// EX: TODO: relearn example syntax
    /// let sam_builder = SamBuilder::default();
    /// let pair_buider = sam_builder.pair_builder();
    /// let pair = pair_builder.cigar1(String::from("1D50M50I50M")).build();
    /// sam_builder.add_pair(pair);
    ///
    pub fn add_pair(&mut self, pair_spec: ReadPairSpec) {
        // TODO: convert cigars into htslib cigars
        assert!(bases1 == "" || quals1 == "" || bases1.len() == quals1.len(), "bases1 and quals1 were different lengths.");
        assert!(bases2 == "" || quals2 == "" || bases2.len() == quals2.len(), "bases2 and quals2 were different lengths.");

        let cigar1 = CigarString::try_from(pair_spec.cigar1.into_bytes().as_ref()).expect("Malformed cigar1");
        let cigar2 = CigarString::try_from(pair_spec.cigar2.into_bytes().as_ref()).expect("Malformed cigar2");
        assert!(unmapped1 || bases1 == "" || bases1.len() == cigar1.clone().into_view(0).end_pos() as usize, "bases1 doesn't agree with cigar on length");
        assert!(unmapped1 || bases2 == "" || bases2.len() == cigar2.clone().into_view(0).end_pos() as usize, "bases2 doesn't agree with cigar on length");
        assert!(unmapped1 || quals1 == "" || quals1.len() == cigar1.clone().into_view(0).end_pos() as usize, "quals1 doesn't agree with cigar on length");
        assert!(unmapped2 || quals2 == "" || quals2.len() == cigar2.clone().into_view(0).end_pos() as usize, "quals2 doesn't agree with cigar on length");


        let mut r1 = Record::new();
        r1.set(
            pair_spec.name.into_bytes().as_ref(),
            if !pair_spec.unmapped1 { Some(cigar1) } else { None },
            pair_spec.bases1.into_bytes().as_ref(),
            pair_spec.quals1.into_bytes().as_ref(),
        );
        r1.set_header(Rc::new(HeaderView::from_header(&self.header)));
        r1.set_tid(pair_spec.contig);
        r1.set_pos(pair_spec.start1);
        if !pair_spec.unmapped1 { r1.set_mapq(pair_spec.mapq1) }
        // TODO: strand
        match pair_spec.strand1 {
            Strand::Plus => (),
            Strand::Minus => r1.set_revere()

        }
        r1.set_paired();
        r1.set_first_in_template();
        if pair_spec.unmapped1 {
            r1.set_unmapped();
        }

        // TODO: set attrs
        // TODO: set readgroup
        // TODO: set mate stuff
        // Add a set mate info method to rust htslib
        // if pair_spec.unmapped2 {
        //     r1.set_mate_unmapped();
        // }


    }
}

impl Default for SamBuilder {
    fn default() -> Self {
        SamBuilder::new(
            DEFAULT_READ_LENGTH,
            DEFAULT_BASE_QUALITY,
            SampleName::default().0,
            Some(ReadGroupId::default().0),
            ReadSortOrder::Unsorted,
            SequenceMetadataStore::default().into(),
            Some(DEFAULT_SEED)
        )
    }
}


/// A specification for a read pair. See ReadPairSpecBuilder or pair_builder
#[derive(Builder)]
pub struct ReadPairSpec {
    name: String,
    bases1: String,
    bases2: String,
    quals1: String,
    quals2: String,
    contig: i32,
    start1: i64,
    start2: i64,
    unmapped1: bool,
    unmapped2: bool,
    cigar1: String,
    cigar2: String,
    mapq1: u8,
    mapq2: u8,
    strand1: Strand,
    strand2: Strand,
    attrs: HashMap<String, String>
}



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
