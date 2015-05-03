# ndarray-blas-level2-complex

# Warning: This package is a placeholder and is not yet tested.

[![Build Status](https://travis-ci.org/scijs/ndarray-blas-level2-complex.svg?branch=master)](https://travis-ci.org/scijs/ndarray-blas-level2-complex) [![npm version](https://badge.fury.io/js/ndarray-blas-level2-complex.svg)](http://badge.fury.io/js/ndarray-blas-level2-complex)

BLAS Level 2 operations for complex [ndarrays](https://github.com/scijs/ndarray)


## Usage

This library implements the basic matrix-vector operations of the Level 2 Basic Linear Algebra Subprograms (BLAS).

Note: It's possible to accomplish the lower triangular functions with the upper triangular version plus flipping and unflipping dimensions, but that's a little convoluted. Instead, the lower triangular versions are suffixed with \_lower just to keep it really simple.

### `gemv( alpha_r, alpha_i, A_r, A_i, x_r, x_i, beta_r, beta_i, y_r, y_i )`
Calculate `y <- alpha*A*x + beta*y`

### `trmv( A_r, A_i, x_r, x_i, isLower )`
Calculate `x <- A*x` for the upper triangular matrix A. Data below the diagonal is ignored. If `isLower` is true, uses the lower triangular portion of A instead.

### `trsv( A_r, A_i, x_r, x_i, isLower )`
Calculate `x <- A^-1 x` for the upper triangular matrix A. Data below the diagonal is ignored.  If `isLower` is true, uses the lower triangular portion of A instead.

## Credits
(c) 2015 Ricky Reusser. MIT License
