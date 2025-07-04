//! Implementation of BFV https://eprint.iacr.org/2012/144.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

pub mod encoder;

pub use encoder::Encoder;
