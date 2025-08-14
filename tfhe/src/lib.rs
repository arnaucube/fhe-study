//! Implementation of TFHE https://eprint.iacr.org/2018/421.pdf
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

pub mod tggsw;
pub mod tglwe;
pub mod tgsw;
pub mod tlev;
pub mod tlwe;

pub(crate) const ERR_SIGMA: f64 = 3.2;
