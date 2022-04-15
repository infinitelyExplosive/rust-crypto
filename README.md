# rust-crypto
Implementing the algorithms and attacks described in _Twenty Years of Attacks on the RSA Cryptosystem_ (Boneh 1999) in Rust. 

Requires [rug](https://crates.io/crates/rug)

## Features
Code is split across 3 files, `main.rs`, which contains code to test algorithms and attacks, `cryptlib.rs`, which contains algorithms used by the attacks, and `cryptlib_bv.rs`, which contains bivariate versions of the algorithms.

### cryptlib
* Univariate polynomial operations
* Determinant
* Resultant
* Euclidean algorithm
* Polynomial euclidean algorithm mod n
* Chinese remainder theorem
* Quadratic equation solver mod n (partial implementation, does not cover all cases)
* Coppersmith's method (Howgrave-Graham simplification)
* LLL algorithm
* Newton's method for approximating zeros

### cryptlib_bv
* Bivariate polynomial operations
* Coppersmith's bivariate method (Coron simplification)
