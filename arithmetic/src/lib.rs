#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

mod naive_ntt; // TODO rm
pub mod ntt;
pub mod ring;
pub mod ringq;
pub mod zq;

pub use ntt::NTT;
pub use ring::R;
pub use ringq::Rq;
pub use zq::Zq;
