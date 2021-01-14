// #![warn(missing_docs)]
use rust_htslib::bam::record::Aux;
/// Wrapper type for a ReadGroupID
#[derive(Debug)]
pub struct ReadGroupId(pub String);

impl Default for ReadGroupId {
    fn default() -> Self {
        Self(String::from("A"))
    }
}

/// Wrapper type for a Sample Name
#[derive(Debug)]
pub struct SampleName(pub String);

impl Default for SampleName {
    fn default() -> Self {
        Self(String::from("Sample1"))
    }
}

/// Wrapper for Aux that owns its own data
#[derive(Clone, Debug)]
pub enum AuxType {
    Integer(i64),
    String(String),
    Float(f64),
    Char(u8),
}

impl<'a> From<&'a AuxType> for Aux<'a> {
    /// Convert an &AuxType to &Aux
    fn from(aux_type: &'a AuxType) -> Self {
        match aux_type {
            AuxType::Integer(x) => Aux::Integer(*x),
            AuxType::String(x) => Aux::String((*x).as_bytes()),
            AuxType::Float(x) => Aux::Float(*x),
            AuxType::Char(x) => Aux::Char(*x),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Strand {
    Plus,
    Minus,
}
