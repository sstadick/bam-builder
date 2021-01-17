//! Simple wrapper over [`rust_htslib::bam`] for building collections of BAM records for testing.
//!
//! Provides an easy to use builder pattern for creating both paired end and single reads.
//! Additionally this library provides nice wrappers and extensions on [rust_htslib] and
//! intends to grow them over time.
//!
//! - [`bam_order`] allows for specifying different sort orders for your BAM collection.
//! - [`sequence_dict`] is a SequenceDictionary like object used for generating [`rust_htslib::bam::Header`]s.
//! - [`wrappers`] provides some convenience types for awkward types like [`rust_htslib::bam::record::Aux`].
//!
//! # Example
//!
//! ```
//! use bam_builder::{bam_order::BamSortOrder, BamBuilder};
//!
//! let mut builder = BamBuilder::new(
//!     100,                    // default read_length
//!     30,                     // default base_quality
//!     "Pair".to_owned(),      // name of samples
//!     None,                   // optional read group id
//!     BamSortOrder::Unsorted, // how to sort reads when `.sort` is called
//!     None,                   // optional sequence_dict
//!     None,                   // seed used for generating random bases
//! );
//!
//! // Create a builder for read pair spec
//!  let records = builder
//!      .pair_builder()
//!      .contig(0)
//!      .start1(0)
//!      .start2(200)
//!      .unmapped1(false)
//!      .unmapped2(false)
//!      .bases1("A".repeat(100))
//!      .bases2("C".repeat(100))
//!      .build()
//!      .unwrap();
//!
//! // Inflate the underlying BAM bam records and add to builder
//! builder.add_pair(records);
//!
//! // Write to temp file
//! let tmp = builder.to_tmp().unwrap();
//! ```
#![warn(missing_docs)]
pub mod bam_order;
pub mod sequence_dict;
pub mod wrappers;
use bam_order::BamSortOrder;
use derive_builder::*;
use rand::prelude::*;
use rust_htslib::bam::{
    header::{Header, HeaderRecord},
    record::{CigarString, Record},
    HeaderView,
};
use rust_htslib::{
    bam::{self, record::Aux},
    errors::Error,
};
use sequence_dict::SequenceDict;
use std::collections::HashMap;
use std::rc::Rc;
use std::{convert::TryFrom, path::Path};
use tempfile::NamedTempFile;
use wrappers::{AuxType, ReadGroupId, SampleName, Strand};

/// Default seed used in random base selection.
const DEFAULT_SEED: usize = 42;
/// Default bases to select from in random base selection.
const DEFAULT_BASES: [char; 4] = ['A', 'C', 'T', 'G'];

/// Builder for BAM records.
///
/// Sets up defaults for reads that are added as well as the BAM header information. Use [`BamBuilder::new`] to generate
/// a [`BamBuilder`] with appropriately populated fields.
///
/// The pattern of interaction with [`BamBuilder`] is to use the [`BamBuilder::frag_builder`] and [`BamBuilder::pair_builder`]
/// methods to create a respective [`ReadFragSpecBuilder`] and [`ReadPairSpecBuilder`]. The `*SpecBuilder`
/// is creates a record with its `.build()` method, and the result can be added to the builder with [`BamBuilder::add_pair`]
/// or [`BamBuilder::add_frag`].
// See: https://github.com/fulcrumgenomics/fgbio/blob/57988d615fa6188e0028faedfdf65ed13cf645e7/src/main/scala/com/fulcrumgenomics/testing/SamBuilder.scala
#[derive(Debug)]
pub struct BamBuilder {
    /// Length of reads to generate
    pub read_length: usize,
    /// Default base quality to use across all qualities
    pub base_quality: u8,
    /// Name of the sample
    pub sample_name: SampleName,
    /// The read group to include in tags
    pub read_group_id: ReadGroupId,
    /// The sort order of the generated reads
    pub sort_order: BamSortOrder,
    /// The BAM header
    pub header: Header,
    /// Random number generator
    pub rng: StdRng,
    /// The collection of records being accumulated
    pub records: Vec<Record>,
    /// The default string o qualities to use
    pub default_quals: String,
    /// Counter for helping to increment the name
    pub counter: usize,
}

impl BamBuilder {
    /// Create a [`BamBuilder`] that fills in defaults based on the passed in options.
    ///
    /// # Example
    /// ```
    /// use bam_builder::{bam_order::BamSortOrder, BamBuilder};
    ///
    /// let mut builder = BamBuilder::new(
    ///     100,                    // default read_length
    ///     30,                     // default base_quality
    ///     "Pair".to_owned(),      // name of samples
    ///     None,                   // optional read group id
    ///     BamSortOrder::Unsorted, // how to sort reads when `.sort` is called
    ///     None,                   // optional sequence_dict
    ///     None,                   // seed used for generating random bases
    /// );
    /// assert_eq!(builder.read_group_id.0, String::from("A"));
    /// assert_eq!(builder.sample_name.0, String::from("Pair"));
    /// assert_eq!(builder.records.len(), 0);
    /// ```
    pub fn new(
        read_length: usize,
        base_quality: u8,
        sample_name: String,
        read_group_id: Option<String>,
        sort: BamSortOrder,
        pseudo_sd: Option<SequenceDict>,
        seed: Option<usize>,
    ) -> BamBuilder {
        let pseudo_sd = pseudo_sd.unwrap_or(SequenceDict::default());
        let seed = seed.unwrap_or(DEFAULT_SEED);
        let rng = StdRng::seed_from_u64(seed as u64);
        let read_group_id = match read_group_id {
            Some(rgi) => ReadGroupId(rgi),
            None => ReadGroupId::default(),
        };
        let default_quals = std::iter::repeat(char::from(33 + base_quality))
            .take(read_length)
            .collect::<String>();

        let mut header: Header = pseudo_sd.into();

        let mut header_record = HeaderRecord::new("RG".as_bytes());
        header_record.push_tag("ID".as_bytes(), &read_group_id.0);
        header.push_record(&header_record);

        BamBuilder {
            read_length,
            base_quality,
            sample_name: SampleName(sample_name),
            read_group_id,
            sort_order: sort,
            header,
            rng,
            records: vec![],
            default_quals,
            counter: 0,
        }
    }

    /// Generate a sequential name, 16 chars wide padded with 0's.
    fn next_name(&mut self) -> String {
        let name = format!("{:0>16}", self.counter);
        self.counter += 1;
        name
    }

    /// Generate a random sequence of bases of length readLength.
    fn random_bases(&mut self) -> String {
        (0..self.read_length)
            .map(|_| DEFAULT_BASES[self.rng.gen_range(0, DEFAULT_BASES.len())])
            .collect::<String>()
    }

    /// Create a ReadFragSpecBuilder with defaults filled in based on BamBuilder.
    pub fn frag_builder(&mut self) -> ReadFragSpecBuilder {
        ReadFragSpecBuilder::default()
            .name(self.next_name())
            .bases(self.random_bases())
            .quals(self.default_quals.clone())
            .contig(-1)
            .start(-1)
            .unmapped(true)
            .cigar(format!("{}M", &self.read_length))
            .mapq(60)
            .strand(Strand::Plus)
            .attrs(HashMap::new())
            .to_owned()
    }

    /// Add a single fragment to the builder.
    pub fn add_frag(&mut self, frag_spec: ReadFragSpec) {
        assert!(
            frag_spec.bases == ""
                || frag_spec.quals == ""
                || frag_spec.bases.len() == frag_spec.quals.len(),
            "bases and quals were different lengths."
        );

        let cigar = CigarString::try_from(frag_spec.cigar.as_str()).expect("Malformed cigar");
        assert!(
            frag_spec.unmapped
                || frag_spec.bases == ""
                || frag_spec.bases.len() == cigar.clone().into_view(0).end_pos() as usize,
            "bases doesn't agree with cigar on length"
        );
        assert!(
            frag_spec.unmapped
                || frag_spec.quals == ""
                || frag_spec.quals.len() == cigar.clone().into_view(0).end_pos() as usize,
            "quals doesn't agree with cigar on length"
        );

        let mut r = Record::new();
        let cigar = CigarString::try_from(frag_spec.cigar.as_str()).unwrap();
        r.set(
            frag_spec.name.into_bytes().as_ref(),
            if !frag_spec.unmapped {
                Some(&cigar)
            } else {
                None
            },
            frag_spec.bases.into_bytes().as_ref(),
            frag_spec.quals.into_bytes().as_ref(),
        );
        r.set_header(Rc::new(HeaderView::from_header(&self.header)));
        r.set_tid(frag_spec.contig);
        r.set_pos(frag_spec.start);
        if !frag_spec.unmapped {
            r.set_mapq(frag_spec.mapq)
        }
        match frag_spec.strand {
            Strand::Plus => (),
            Strand::Minus => r.set_reverse(),
        }
        if frag_spec.unmapped {
            r.set_unmapped();
        }
        r.push_aux(
            "RG".as_bytes(),
            &Aux::String(&self.read_group_id.0.as_bytes()),
        );
        for (key, value) in frag_spec.attrs.iter() {
            r.push_aux(key.as_bytes(), &value.into());
        }
        self.records.push(r);
    }

    /// Create a ReadPairSpecBuilder with defaults filled in based on the BamBuilder.
    pub fn pair_builder(&mut self) -> ReadPairSpecBuilder {
        ReadPairSpecBuilder::default()
            .name(self.next_name())
            .bases1(self.random_bases())
            .bases2(self.random_bases())
            .quals1(self.default_quals.clone())
            .quals2(self.default_quals.clone())
            .contig(-1)
            .start1(-1) // corresponds to unset
            .start2(-1)
            .unmapped1(true)
            .unmapped2(true)
            .cigar1(format!("{}M", &self.read_length))
            .cigar2(format!("{}M", &self.read_length))
            .mapq1(60)
            .mapq2(60)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .attrs(HashMap::new())
            .to_owned()
    }

    /// Add a read pair to the builder.
    pub fn add_pair(&mut self, pair_spec: ReadPairSpec) {
        assert!(
            pair_spec.bases1 == ""
                || pair_spec.quals1 == ""
                || pair_spec.bases1.len() == pair_spec.quals1.len(),
            "bases1 and quals1 were different lengths."
        );
        assert!(
            pair_spec.bases2 == ""
                || pair_spec.quals2 == ""
                || pair_spec.bases2.len() == pair_spec.quals2.len(),
            "bases2 and quals2 were different lengths."
        );

        let cigar1 = CigarString::try_from(pair_spec.cigar1.as_str()).expect("Malformed cigar1");
        let cigar2 = CigarString::try_from(pair_spec.cigar2.as_str()).expect("Malformed cigar2");
        assert!(
            pair_spec.unmapped1
                || pair_spec.bases1 == ""
                || pair_spec.bases1.len() == cigar1.clone().into_view(0).end_pos() as usize,
            "bases1 doesn't agree with cigar on length"
        );
        assert!(
            pair_spec.unmapped1
                || pair_spec.bases2 == ""
                || pair_spec.bases2.len() == cigar2.clone().into_view(0).end_pos() as usize,
            "bases2 doesn't agree with cigar on length"
        );
        assert!(
            pair_spec.unmapped1
                || pair_spec.quals1 == ""
                || pair_spec.quals1.len() == cigar1.clone().into_view(0).end_pos() as usize,
            "quals1 doesn't agree with cigar on length"
        );
        assert!(
            pair_spec.unmapped2
                || pair_spec.quals2 == ""
                || pair_spec.quals2.len() == cigar2.clone().into_view(0).end_pos() as usize,
            "quals2 doesn't agree with cigar on length"
        );

        let mut r1 = Record::new();
        let cigar = CigarString::try_from(pair_spec.cigar1.as_str()).unwrap();
        r1.set(
            pair_spec.name.into_bytes().as_ref(),
            if !pair_spec.unmapped1 {
                Some(&cigar)
            } else {
                None
            },
            pair_spec.bases1.into_bytes().as_ref(),
            pair_spec.quals1.into_bytes().as_ref(),
        );
        r1.set_header(Rc::new(HeaderView::from_header(&self.header)));
        r1.set_tid(pair_spec.contig);
        r1.set_pos(pair_spec.start1);
        if !pair_spec.unmapped1 {
            r1.set_mapq(pair_spec.mapq1)
        }
        match pair_spec.strand1 {
            Strand::Plus => (),
            Strand::Minus => r1.set_reverse(),
        }
        r1.set_paired();
        r1.set_first_in_template();
        if pair_spec.unmapped1 {
            r1.set_unmapped();
        }

        let mut r2 = Record::new();
        let cigar = CigarString::try_from(pair_spec.cigar2.as_str()).unwrap();
        r2.set(
            r1.qname(),
            if !pair_spec.unmapped2 {
                Some(&cigar)
            } else {
                None
            },
            pair_spec.bases2.into_bytes().as_ref(),
            pair_spec.quals2.into_bytes().as_ref(),
        );
        r2.set_header(Rc::new(HeaderView::from_header(&self.header)));
        r2.set_tid(pair_spec.contig);
        r2.set_pos(pair_spec.start2);
        if !pair_spec.unmapped2 {
            r2.set_mapq(pair_spec.mapq2)
        }
        match pair_spec.strand2 {
            Strand::Plus => (),
            Strand::Minus => r2.set_reverse(),
        }
        r2.set_paired();
        r2.set_first_in_template();
        if pair_spec.unmapped2 {
            r2.set_unmapped();
        }
        r1.push_aux(
            "RG".as_bytes(),
            &Aux::String(&self.read_group_id.0.as_bytes()),
        );
        r2.push_aux(
            "RG".as_bytes(),
            &Aux::String(&self.read_group_id.0.as_bytes()),
        );

        for (key, value) in pair_spec.attrs.iter() {
            r1.push_aux(key.as_bytes(), &value.into());
            r2.push_aux(key.as_bytes(), &value.into());
        }
        BamBuilder::set_mate_info(&mut r1, &mut r2, true);
        self.records.push(r1);
        self.records.push(r2);
    }

    /// Write the mate info for two BAM Records.
    fn set_mate_info(rec1: &mut Record, rec2: &mut Record, set_mate_cigar: bool) {
        // See SamPairUtil in htsjdk for original impl
        // If neither read is unmapped just set their mate info
        if !rec1.is_unmapped() && !rec2.is_unmapped() {
            rec1.set_mtid(rec2.tid());
            rec1.set_mpos(rec2.pos());
            if rec2.is_reverse() {
                rec1.set_mate_reverse()
            }
            rec1.push_aux(b"MQ", &Aux::Char(rec2.mapq()));

            rec2.set_mtid(rec1.tid());
            rec2.set_mpos(rec1.pos());
            if rec1.is_reverse() {
                rec2.set_mate_reverse()
            }
            rec2.push_aux(b"MQ", &Aux::Char(rec1.mapq()));

            if set_mate_cigar {
                rec1.push_aux(b"MC", &Aux::String(rec2.cigar().to_string().as_bytes()));
                rec2.push_aux(b"MC", &Aux::String(rec1.cigar().to_string().as_bytes()));
            } // leave empty otherwise
        } else if rec1.is_unmapped() && rec2.is_unmapped() {
            // both unmapped
            rec1.set_tid(-1);
            rec1.set_pos(-1); // corresponds to unset
            rec1.set_mtid(-1);
            rec1.set_mpos(-1);
            rec1.set_mate_unmapped();
            rec1.set_insert_size(0);
            if rec2.is_reverse() {
                rec1.set_mate_reverse()
            }

            rec2.set_tid(-1);
            rec2.set_pos(-1); // corresponds to unset
            rec2.set_mtid(-1);
            rec2.set_mpos(-1);
            rec2.set_mate_unmapped();
            rec2.set_insert_size(0);
            if rec1.is_reverse() {
                rec2.set_mate_reverse()
            }
        } else if rec1.is_unmapped() {
            // only rec2 is mapped
            rec1.set_tid(-1);
            rec1.set_pos(-1); // corresponds to unset

            rec2.set_mtid(rec1.tid());
            rec2.set_mpos(rec1.pos());
            rec2.set_mate_unmapped();
            rec2.set_insert_size(0);
            if rec1.is_reverse() {
                rec2.set_mate_reverse()
            }

            rec1.set_mtid(rec2.tid());
            rec1.set_mpos(rec2.pos());
            rec1.set_insert_size(0);
            rec1.push_aux(b"MQ", &Aux::Char(rec2.mapq()));
            if set_mate_cigar {
                rec1.push_aux(b"MC", &Aux::String(rec2.cigar().to_string().as_bytes()))
            }
            if rec2.is_reverse() {
                rec1.set_mate_reverse()
            }
        } else {
            // only rec1 is mapped
            rec2.set_tid(-1);
            rec2.set_pos(-1); // corresponds to unset

            rec1.set_mtid(rec2.tid());
            rec1.set_mpos(rec2.pos());
            rec1.set_mate_unmapped();
            rec1.set_insert_size(0);
            if rec2.is_reverse() {
                rec1.set_mate_reverse()
            }

            rec2.set_mtid(rec1.tid());
            rec2.set_mpos(rec1.pos());
            rec2.set_insert_size(0);
            rec2.push_aux(b"MQ", &Aux::Char(rec1.mapq()));
            if set_mate_cigar {
                rec2.push_aux(b"MC", &Aux::String(rec1.cigar().to_string().as_bytes()))
            }
            if rec1.is_reverse() {
                rec2.set_mate_reverse()
            }
        }

        let insert_size = BamBuilder::compute_insert_size(rec1, rec2);
        rec1.set_insert_size(insert_size);
        rec2.set_insert_size(insert_size);
    }

    /// Computes the TLEN field for two records.
    fn compute_insert_size(rec1: &Record, rec2: &Record) -> i64 {
        if rec1.is_unmapped() || rec2.is_unmapped() {
            return 0;
        }
        if rec1.tid() != rec2.tid() {
            return 0;
        }

        let rec1_5prime_pos = if rec1.is_reverse() {
            rec1.cigar().end_pos()
        } else {
            rec1.pos()
        };
        let rec2_5prime_pos = if rec2.is_reverse() {
            rec2.cigar().end_pos()
        } else {
            rec2.pos()
        };

        rec2_5prime_pos - rec1_5prime_pos //+ adjustment
    }

    /// Writes BAM records to specified path.
    ///
    /// Note that [`BamBuilder::sort`] should be called ahead of writing to ensure records
    /// are sorted.
    pub fn to_path(&self, path: &Path) -> Result<(), Error> {
        let mut writer = bam::Writer::from_path(path, &self.header, bam::Format::BAM)?;
        for record in self.records.iter() {
            writer.write(record)?;
        }

        Ok(())
    }

    /// Write records to a tempfile. Tempfile will be deleted when `NamedTempFile` goes out of scope.
    ///
    /// Note that [`BamBuilder::sort`] should be called ahead of writing to ensure records
    /// are sorted.
    pub fn to_tmp(&self) -> Result<NamedTempFile, Error> {
        let tempfile = NamedTempFile::new().unwrap();
        self.to_path(tempfile.path())?;
        Ok(tempfile)
    }

    /// Sort the records added thus far by whichever [`BamSortOrder`] was specified.
    pub fn sort(&mut self) {
        self.sort_order.sort(&mut self.records);
    }
}

/// A specification for a read pair. See [`ReadPairSpecBuilder`] or [`BamBuilder::pair_builder`].
/// Mate tags and info will be set by default
#[derive(Builder, Debug)]
pub struct ReadPairSpec {
    /// Name of the reads, will be filled in by [`BamBuilder::next_name`] by default.
    name: String,
    /// Sequence bases, will be filled in by [`BamBuilder::random_bases`] by default.
    bases1: String,
    /// Sequence bases, will be filled in by [`BamBuilder::random_bases`] by default.
    bases2: String,
    /// Quality string, will be based on the [`BamBuilder::base_quality`] by default.
    quals1: String,
    /// Quality string, will be based on the [`BamBuilder::base_quality`] by default.
    quals2: String,
    /// Reference contig if mapped, defaults to unmapped.
    contig: i32,
    /// 0-based location on reference contig if mapped. defaults to unmapped.
    start1: i64,
    /// 0-based location on reference contig if mapped. defaults to unmapped.
    start2: i64,
    /// `true` if unmapped, defaults to `true`.
    unmapped1: bool,
    /// `true` if unmapped, defaults to `true`.
    unmapped2: bool,
    /// Alignment, defaults to all `M` if mapped or `*` if unmapped.
    cigar1: String,
    /// Alignment, defaults to all `M` if mapped or `*` if unmapped.
    cigar2: String,
    /// map quality of the read, defaults to 60 if mapped, `*` if unmapped.
    mapq1: u8,
    /// map quality of the read, defaults to 60 if mapped, `*` if unmapped.
    mapq2: u8,
    /// Strand the read maps to, defaults to [`wrappers::Strand::Plus`].
    strand1: Strand,
    /// Strand the read maps to, defaults to [`wrappers::Strand::Minus`].
    strand2: Strand,
    /// Tags for the reads, mate tags will be filled in by default.
    attrs: HashMap<String, AuxType>,
}

/// A specification for a single read. See [`ReadFragSpecBuilder`] or [`BamBuilder::frag_builder`].
#[derive(Builder, Debug)]
pub struct ReadFragSpec {
    /// Name of the read, will be filled in by [`BamBuilder::next_name`] by default.
    name: String,
    /// Sequence bases, will be filled in by [`BamBuilder::random_bases`] by default.
    bases: String,
    /// Quality string, will be based on the [`BamBuilder::base_quality`] by default.
    quals: String,
    /// Reference contig if mapped, defaults to unmapped.
    contig: i32,
    /// 0-based location on reference contig if mapped. defaults to unmapped.
    start: i64,
    /// `true` if unmapped, defaults to `true`.
    unmapped: bool,
    /// Alignment, defaults to all `M` if mapped or `*` if unmapped.
    cigar: String,
    /// map quality of the read, defaults to 60 if mapped, `*` if unmapped.
    mapq: u8,
    /// Strand the read maps to, defaults to [`wrappers::Strand::Plus`].
    strand: Strand,
    /// Tags for the read, defaults to none.
    attrs: HashMap<String, AuxType>,
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
