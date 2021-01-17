//! Helper wrappers over awkward [`rust_htslib`] types.
#![warn(missing_docs)]
use rust_htslib::bam::record::Aux;
/// Wrapper type for a ReadGroupID.
#[derive(Debug)]
pub struct ReadGroupId(pub String);

impl Default for ReadGroupId {
    /// Default [`ReadGroupId`] is "A".
    fn default() -> Self {
        Self(String::from("A"))
    }
}

/// Wrapper type for a Sample Name
#[derive(Debug)]
pub struct SampleName(pub String);

impl Default for SampleName {
    /// Default [`SampleName`] is "Sample1".
    fn default() -> Self {
        Self(String::from("Sample1"))
    }
}

/// Wrapper for Aux that owns its own data. [`Aux::from`] is implemented in order
/// to facilitate conversion between the wrapper and the underlying [`Aux`] type.
#[derive(Clone, Debug)]
pub enum AuxType {
    /// Corresponds to [Aux::Integer].
    Integer(i64),
    /// Corresponds to [Aux::String], but as an owned String.
    String(String),
    /// Corresponds to [Aux::Float].
    Float(f64),
    /// Corresponds to [Aux::Char].
    Char(u8),
}

impl<'a> From<&'a AuxType> for Aux<'a> {
    /// Convert a ref [`AuxType`] to ref [`Aux`].
    fn from(aux_type: &'a AuxType) -> Self {
        match aux_type {
            AuxType::Integer(x) => Aux::Integer(*x),
            AuxType::String(x) => Aux::String((*x).as_bytes()),
            AuxType::Float(x) => Aux::Float(*x),
            AuxType::Char(x) => Aux::Char(*x),
        }
    }
}

/// Wrapper around strandedness of a read.
#[derive(Clone, Debug)]
pub enum Strand {
    /// Forward strand.
    Plus,
    /// Reverse strand.
    Minus,
}
