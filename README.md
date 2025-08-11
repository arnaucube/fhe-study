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
and the line `type S = TWLE<K>` to use `CKKS<Q, N>` or `BFV<Q, N, T>`.

```rust
const T: u64 = 128; // msg space (msg modulus)
type M = Rq<T, 1>; // msg space
type S = TLWE<256>;

let mut rng = rand::thread_rng();
let msg_dist = Uniform::new(0_u64, T);

let (sk, pk) = S::new_key(&mut rng)?;

// get two random msgs in Z_t
let m1 = M::rand_u64(&mut rng, msg_dist)?;
let m2 = M::rand_u64(&mut rng, msg_dist)?;
let m3 = M::rand_u64(&mut rng, msg_dist)?;

// encode the msgs into the plaintext space
let p1 = S::encode::<T>(&m1); // plaintext
let p2 = S::encode::<T>(&m2); // plaintext
let c3_const: Tn<1> = Tn(array::from_fn(|i| T64(m3.coeffs()[i].0))); // encode it as constant

let c1 = S::encrypt(&mut rng, &pk, &p1)?;
let c2 = S::encrypt(&mut rng, &pk, &p2)?;

// now we can do encrypted operations (notice that we do them using simple
// operation notation by rust's operator overloading):
let c_12 = c1 + c2;
let c4 = c_12 * c3_const;

// decrypt & decode
let p4_recovered = c4.decrypt(&sk);
let m4 = S::decode::<T>(&p4_recovered);

// m4 is equal to (m1+m2)*m3
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
