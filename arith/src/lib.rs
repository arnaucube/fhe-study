#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(dead_code)] // TMP

pub mod complex;
pub mod matrix;
pub mod torus;
pub mod zq;

pub mod ring;
pub mod ring_n;
pub mod ring_nq;
pub mod ring_torus;
pub mod tuple_ring;

// mod naive_ntt; // note: for dev only
pub mod ntt;
pub mod ntt_u62;
pub mod ntt_u64;

// expose objects
pub use complex::C;
pub use matrix::Matrix;
pub use torus::T64;
pub use zq::Zq;

pub use ring::{Ring, RingParam};
pub use ring_n::R;
pub use ring_nq::Rq;
pub use ring_torus::Tn;
pub use tuple_ring::TR;

pub use ntt::NTT;
