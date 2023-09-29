//! Wrapper over [`rust_htslib::bam::header`], see [`SequenceDict`].
//!
//! The intention is to turn this into something like Picards's SequenceDictionary at
//! some future point.
#![warn(missing_docs)]
use rust_htslib::bam::header::{Header, HeaderRecord};

/// Hold relevant info about the sequences, translates into a header.
///
/// # Example
/// ```
/// use bam_builder::sequence_dict::{SequenceDict, SequenceData};
/// let sequences = (1..22)
/// .map(|chr| chr.to_string())
/// .chain(vec!["X", "Y", "M"].into_iter().map(|chr| chr.to_string()))
/// .map(|chr| SequenceData {
///     name: format!("chr{}", chr),
///     length: 2_000_000,
/// })
/// .collect();
/// let sd = SequenceDict { sequences };
/// assert_eq!(sd.sequences[0].name, String::from("chr1"));
/// ```
#[derive(Debug)]
pub struct SequenceDict {
    /// Ordered list of [`SequenceData`].
    pub sequences: Vec<SequenceData>,
}

/// Holds relevant information about reference sequences
// Note: At a future point spin this out into its own lib and support Sequence Dict
#[derive(Debug)]
pub struct SequenceData {
    /// Name of the sequence.
    pub name: String,
    /// Length of the sequence.
    pub length: usize,
}

impl From<SequenceData> for HeaderRecord<'_> {
    fn from(meta: SequenceData) -> Self {
        let mut rec = HeaderRecord::new("SQ".as_bytes());
        rec.push_tag("SN".as_bytes(), &meta.name);
        rec.push_tag("LN".as_bytes(), meta.length);
        rec
    }
}

impl Default for SequenceDict {
    fn default() -> Self {
        let sequences = (1..22)
            .map(|chr| chr.to_string())
            .chain(vec!["X", "Y", "M"].into_iter().map(|chr| chr.to_string()))
            .map(|chr| SequenceData {
                name: format!("chr{}", chr),
                length: 2_000_000,
            })
            .collect();
        SequenceDict { sequences }
    }
}

impl From<SequenceDict> for Header {
    fn from(meta: SequenceDict) -> Header {
        let mut header = Header::new();
        for seq in meta.sequences.into_iter() {
            header.push_record(&seq.into());
        }
        header
    }
}
