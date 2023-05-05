#include <Rcpp.h>
using namespace Rcpp;

#include <cstdlib>
#include <iostream>
#include "miniball.hpp"

// Rcpp miniball function, adapted from code: http://people.inf.ethz.ch/gaertner/subdir/software/miniball/miniball_example.cpp

// [[Rcpp::export]]
double miniball_rsqr(NumericMatrix xmat, int d, int n) {

  typedef double mytype;            // coordinate type

  // populate 2d array
  // ----------------------------------------------------
  mytype** ap = new mytype*[n];
  for (int i=0; i<n; ++i) {
    mytype* p = new mytype[d];
    for (int j=0; j<d; ++j) {
      p[j] = xmat[(i*d)+j];
    }
    ap[i]=p;
  }

  // iterator types
  typedef mytype* const* PointIterator;
  typedef const mytype* CoordIterator;
  float rsqr;
  
  // instance of miniball
  // d is latent dimension, n is number of points, ap is n x d matrix of points
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> >
    MB;
  MB mb (d, ap, ap+n);

  rsqr = mb.squared_radius();
  
  return rsqr;
}
