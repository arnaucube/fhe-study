#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

pub mod complex;
pub mod matrix;
mod naive_ntt; // note: for dev only
pub mod ntt;
pub mod ring;
pub mod ringq;
pub mod zq;

pub use complex::C;
pub use matrix::Matrix;
pub use ntt::NTT;
pub use ring::R;
pub use ringq::Rq;
pub use zq::Zq;
