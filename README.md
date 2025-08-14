# fhe-study
Implementations from scratch done while studying some FHE papers; do not use in production.

- `arith`: contains $\mathbb{Z}_q$, $R_q=\mathbb{Z}_q[X]/(X^N+1)$, $R=\mathbb{Z}[X]/(X^N+1)$, $\mathbb{T}_q[X]/(X^N +1)$ arithmetic implementations, together with the NTT implementation.
- `gfhe`: (gfhe=generalized-fhe) contains the structs and logic for RLWE, GLWE, GLev, GGSW, RGSW cryptosystems, and modulus switching and key switching methods, which can be used by concrete FHE schemes.
- `bfv`: https://eprint.iacr.org/2012/144.pdf scheme implementation
- `ckks`: https://eprint.iacr.org/2016/421.pdf scheme implementation
- `tfhe`: https://eprint.iacr.org/2018/421.pdf scheme implementation


## Run tests
`cargo test --release`

## Example of usage
> the repo is a work in progress, interfaces will change.

This example shows usage of TFHE, but the idea is that the same interface would
work for using CKKS & BFV, the only thing to be changed would be the parameters
and the usage of `TLWE` by `CKKS` or `BFV`.

```rust
let param = Param {
    err_sigma: crate::ERR_SIGMA,
    ring: RingParam { q: u64::MAX, n: 1 },
    k: 256,
    t: 128, // plaintext modulus
};

let mut rng = rand::thread_rng();
let msg_dist = Uniform::new(0_u64, param.t);

let (sk, pk) = TLWE::new_key(&mut rng, &param)?;

// get three random msgs in Rt
let m1 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
let m2 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;
let m3 = Rq::rand_u64(&mut rng, msg_dist, &param.pt())?;

// encode the msgs into the plaintext space
let p1 = TLWE::encode(&param, &m1); // plaintext
let p2 = TLWE::encode(&param, &m2); // plaintext
let c3_const = TLWE::new_const(&param, &m3); // as constant/public value

// encrypt p1 and m2
let c1 = TLWE::encrypt(&mut rng, &param, &pk, &p1)?;
let c2 = TLWE::encrypt(&mut rng, &param, &pk, &p2)?;

// now we can do encrypted operations (notice that we do them using simple
// operation notation by rust's operator overloading):
let c_12 = c1 + c2;
let c4 = c_12 * c3_const;

// decrypt & decode
let p4_recovered = c4.decrypt(&sk);
let m4 = TLWE::decode(&param, &p4_recovered);

// m4 is equal to (m1+m2)*m3
assert_eq!(((m1 + m2).to_r() * m3.to_r()).to_rq(param.t), m4);
```


## Status of the implementation

- TFHE
	- {TLWE, TGLWE, TLev, TGLev, TGSW, TGGSW} encryption & decryption
	- addition of ciphertexts, addition & multiplication of ciphertext by a plaintext
	- external products of ciphertexts
		- TGSW x TLWE
		- TGGSW x TGLWE
	- {TGSW, TGGSW} CMux gate
	- blind rotation, key switching, mod switching
	- bootstrapping
- CKKS
	- encoding & decoding
	- encryption & decryption
	- addition & substraction of ciphertexts
- BFV
	- encryption & decryption
	- addition & substraction of ciphertexts
	- addition & multiplication of ciphertext by a plaintext
	- multiplication of ciphertexts with relinearization
- GFHE (generalized FHE)
	- {GLWE & GLev} encryption & decryption
	- key switching, mod switching
	- addition & substraction of ciphertexts
	- addition & multiplication of ciphertext by a plaintext
- arith
	- base arithmetic for $\mathbb{Z}_q,~~ R_q=\mathbb{Z}_q[X]/(X^N+1),~~ R=\mathbb{Z}[X]/(X^N+1),~~ \mathbb{T}_q[X]/(X^N +1)$
	- NTT implementation (negacyclic convolution)
