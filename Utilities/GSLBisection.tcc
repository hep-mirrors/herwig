// -*- C++ -*-
//
// GSLBisection.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the GSLBisection class.
//
using namespace Herwig;
using namespace ThePEG;

namespace {

  template<class T> double func(double x , void * p) {
    //Units of the argument and return type
    const typename T::ValType ValUnit = ((struct GSLparam<T> *)p)->function.vUnit();
    const typename T::ArgType ArgUnit = ((struct GSLparam<T> *)p)->function.aUnit();

    const T & f = ((struct GSLparam<T> *)p)->function;
    return f(x * ArgUnit ) / ValUnit;
  }
  
}

namespace Herwig {
using namespace ThePEG;

template <class T> inline typename T::ArgType 
GSLBisection::value(const T & fn, 
		    const typename T::ArgType lower, 
		    const typename T::ArgType upper) const {

  typedef typename T::ArgType ArgType;
  const ArgType ArgUnit = fn.aUnit();
  
  //use own error handler
  gsl_error_handler_t *old_handler = 
    gsl_set_error_handler(& GSLsubstHandler);


  int status(0), iter(0);
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *solver;
  double result(0);
  double x_lo(lower/ArgUnit), x_hi(upper/ArgUnit);

  GSLparam<T> parameters = { fn };
  gsl_function F;
  F.function = & func<T>;
  F.params = &parameters;
     
  solverType = gsl_root_fsolver_brent;
  solver = gsl_root_fsolver_alloc (solverType);

  try{
    gsl_root_fsolver_set (solver, &F, x_lo, x_hi);
  }catch(GSLerror){
    //cerr << "GSLBisection: initial interval does not contain zero\n";
    throw IntervalError();
  }
   
/*    
  printf ("Root finding is using %s method\n", 
	  gsl_root_fsolver_name (solver));
  printf ("%5s [%9s, %9s] %9s %10s\n",
	  "iter", "lower", "upper", "root", 
	  "err");
*/   
  do{
    iter++;
    status = gsl_root_fsolver_iterate (solver);
    result = gsl_root_fsolver_root (solver);
    x_lo = gsl_root_fsolver_x_lower (solver);
    x_hi = gsl_root_fsolver_x_upper (solver);
    status = gsl_root_test_interval (x_lo, x_hi, abserr_, relerr_);
    
/*    
    if (status == GSL_SUCCESS)
      printf ("Converged:\n");

    printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
	    iter, x_lo, x_hi,
	    result, x_hi - x_lo);
*/
  }
  while (status == GSL_CONTINUE && iter < maxPoints_);
    
  gsl_root_fsolver_free (solver);
  //use default GSL error handler again
  gsl_set_error_handler(old_handler);

  //fix units and return
  return result * ArgUnit;
}

}
