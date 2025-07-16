# fhe-study
Implementations from scratch done while studying some FHE papers; do not use in production.

- `arith`: contains $\mathbb{Z}_q$, $R_q=\mathbb{Z}_q[X]/(X^N+1)$ and $R=\mathbb{Z}[X]/(X^N+1)$ arithmetic implementations, together with the NTT implementation.
- `gfhe`: (gfhe=generalized-fhe) contains the structs and logic for RLWE, GLWE, GLev, GGSW, RGSW cryptosystems, and modulus switching and key switching methods, which can be used by concrete FHE schemes.
- `bfv`: https://eprint.iacr.org/2012/144.pdf scheme implementation
- `ckks`: https://eprint.iacr.org/2016/421.pdf scheme implementation

`cargo test --release`
