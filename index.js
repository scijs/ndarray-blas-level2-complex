'use strict';

var blas1 = require('ndarray-blas-level1-complex');

exports.gemv = function(alpha_r, alpha_i, A_r, A_i, x_r, x_i, beta_r, beta_i, y_r, y_i) {
  var dotu = blas1.dotu;
  var Ax, alphaAx_r, alphaAx_i, yBeta_r, yBeta_i, yi_r, yi_i;
  for(var i=A_r.shape[1]-1; i>=0; i--) {
    Ax = dotu( A_r.pick(i,null), A_i.pick(i,null), x_r, x_i);
    alphaAx_r = alpha_r * Ax[0] - alpha_i * Ax[1];
    alphaAx_i = alpha_r * Ax[1] + alpha_i * Ax[0];
    yi_r = y_r.get(i);
    yi_i = y_i.get(i);
    yBeta_r = yi_r * beta_r - yi_i * beta_i;
    yBeta_i = yi_r * beta_i + yi_i * beta_r;

    y_r.set(i, yBeta_r + alphaAx_r);
    y_i.set(i, yBeta_i + alphaAx_i);
  }
  return true;
};

exports.gbmv = function() {
  console.error('GBMV (banded matrix vector multiply) not yet implemented');
};

exports.symv = function() {
  console.error('SYMV (symmetric matrix vector multiply) not yet implemented');
};

exports.sbmv = function() {
  console.error('SBMV (symmetric banded matrix vector multiply) not yet implemented');
};

exports.spmv = function() {
  console.error('SPMV (symmetric packed matrix vector multiply) not yet implemented');
};


// Compute the product of an upper triangular matrix with a vector
exports.trmv = function(A_r, A_i, x_r, x_i, isLower) {
  var dotu = blas1.dotu;
  var n = A_r.shape[1];
  var d;
  if( isLower ) {
    for(var i=n-1; i>=0; i--) {
      d = dotu( A_r.pick(i,null).hi(i+1), A_i.pick(i,null).hi(i+1), x_r.hi(i+1), x_i.hi(i+1) );
      x_r.set(i, d[0] );
      x_i.set(i, d[1] );
    }
  } else {
    for(var i=0; i<n; i++) {
      d = dotu( A_r.pick(i,null).lo(i), A_r.pick(i,null).lo(i), x_r.lo(i), x_i.lo(i) );
      x_r.set(i, d[0] );
      x_i.set(i, d[1] );
    }
  }
  return true;
};

exports.trmv_lower = function(A,x) {
  console.warn('trmv_lower is deprected. Please use the \'isLower\' flag with trmv');
  return exports.trmv(A,x,true);
}

exports.tbmv = function() {
  console.error('TBMV (triangular banded matrix vector multiply) not yet implemented');
};

// Solve Ax=b where A is upper triangular
exports.trsv = function(A_r, A_i, x_r, x_i, isLower) {

  console.error('trsv for triangular matrices not yet implemented');
  return false;

  var i;
  var dotu = blas1.dotu;
  var n = A_r.shape[1];
  var d, x0_r, x0_i, A00_r, A00_i;;
  if( isLower ) {

    x0_r = x_r.get(0);
    x0_i = x_i.get(0);
    A00_r = A_r.get(0,0);
    A00_i = A_i.get(0,0);
    d = A00_r * A00_r + A00_i * A00_i;

    x_r.set( 0, (A00_r * x0_r + A00_i * x0_i) / d );
    x_i.set( 0, (A00_r * x0_i - A00_i * x0_r) / d );

    for(i=1; i<n; i++) {
      d = dotu(A_r.pick(i,null).hi(i), A_i.pick(i,null).hi(i), x_r.hi(i), x_i.hi(i));
      x.set(i, (x.get(i) - d) / A.get(i,i) );
    }
  } else {
    //x.set( n-1, x.get(n-1)/A.get(n-1,n-1) );
    //for(i=n-2; i>=0; i--) {
      //x.set(i, (x.get(i) - dotu(A.pick(i,null).lo(i+1), x.lo(i+1))) / A.get(i,i) );
    //}
  }
  return true;
};

// Solve Ax=b where A is lower triangular
exports.trsv_lower = function(A, x) {
  console.warn('trsv_lower is deprected. Please use the \'isLower\' flag with trsv');
  return exports.trsv(A,x,true);
};

exports.tbsv = function() {
  console.error('TBSV (triangular banded matrix solver) not yet implemented');
};

exports.tpsv = function() {
  console.error('TPSV (triangular packed matrix solver) not yet implemented');
};

exports.ger = function() {
  console.error('GER (rank 1 operation A := alpha*x*y\' + A) not yet implemented');
};

exports.syr = function() {
  console.error('SYR (symmetric rank 1 operation A := alpha*x*y\' + A) not yet implemented');
};

exports.spr = function() {
  console.error('SPR (symmetric packed rank 1 operation A := alpha*x*y\' + A) not yet implemented');
};

exports.syr2 = function() {
  console.error('SYR (symmetric rank 2 operation A := alpha*x*y\' + alpha*y*x\' + A) not yet implemented');
};

exports.spr2 = function() {
  console.error('SPR (symmetric packed rank 2 operation A := alpha*x*y\' + alpha*y*x\' + A) not yet implemented');
};
