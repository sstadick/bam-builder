//! Specification for how to sort BAM records based on [`BamSortOrder`].
#![warn(missing_docs)]
use super::DEFAULT_SEED;
use rand::prelude::*;
use rust_htslib::bam::record::Record;
use std::cmp::Ordering;
/// Sort order options for reads
#[derive(Debug)]
pub enum BamSortOrder {
    /// Sort by `qname`.
    NameSorted,
    /// Sort by `tid` and `pos`.
    CoordSorted,
    /// Leave records in current order.
    Unsorted,
    /// Randomly shuffle records.
    Random,
}

impl BamSortOrder {
    /// Sort `records` based on the [`BamSortOrder`] variant.
    pub fn sort(&self, records: &mut [Record]) {
        match self {
            BamSortOrder::NameSorted => records.sort_by(|a, b| {
                let name_cmp = a.qname().cmp(b.qname());
                match name_cmp {
                    Ordering::Equal => a.flags().cmp(&b.flags()),
                    _ => name_cmp,
                }
            }),
            BamSortOrder::CoordSorted => records.sort_by(|a, b| {
                let ref_cmp = a.tid().cmp(&b.tid());
                match ref_cmp {
                    Ordering::Equal => {
                        let pos_cmp = a.pos().cmp(&b.pos());
                        match pos_cmp {
                            Ordering::Equal => a.flags().cmp(&b.flags()),
                            _ => pos_cmp,
                        }
                    }
                    _ => ref_cmp,
                }
            }),
            BamSortOrder::Unsorted => {} // do nothing
            BamSortOrder::Random => {
                let mut rng = StdRng::seed_from_u64(DEFAULT_SEED as u64);
                records.shuffle(&mut rng);
            }
        }
    }
}

impl Default for BamSortOrder {
    /// Default for [`BamSortOrder`] is [`BamSortOrder::Unsorted`].
    fn default() -> Self {
        Self::Unsorted
    }
}
