// -*- C++ -*-

// math/gsl.tcc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include <cassert>

#include "gsl/gsl_errno.h"

namespace matchbox { namespace gsl {

template<class Function>
double function_wrapper (double x, void * fptr) {
  return (*reinterpret_cast<Function *>(fptr)).operator () (x);
}

/// given a unary function, return an appropriate
/// gsl function object
template<class Function>
gsl_function make_function (Function& f) {
  gsl_function conv;
  conv.function = &function_wrapper<Function>;
  conv.params = &f;
  return conv;
}

template<class Function, unsigned long MaxIterations>
double bisection_root_solver<Function,MaxIterations>::solve (std::pair<double,double> interval, double precision) {

  assert(interval.first < interval.second);
  
  gsl_function F = make_function(f);
  
  gsl_root_fsolver_set (s, &F, interval.first, interval.second);
  
  unsigned long iterations = 0;
  double sol;
  int status;
  
  do {
    ++iterations;
    status = gsl_root_fsolver_iterate (s);
    sol = gsl_root_fsolver_root (s);
    interval.first = gsl_root_fsolver_x_lower(s);
    interval.second = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval (interval.first,interval.second,0,precision);
  } while (status == GSL_CONTINUE && iterations < MaxIterations);
  
  return sol;
  
}

  }}
