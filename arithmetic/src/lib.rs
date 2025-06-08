#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

mod naive; // TODO rm
pub mod ntt;
pub mod ring;
pub mod zq;

pub use ntt::NTT;
pub use ring::PR;
pub use zq::Zq;
