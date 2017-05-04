// -*- C++ -*-
//
// GSLIntegrator.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the GSLIntegrator class.
//
using namespace Herwig;
using namespace ThePEG;

namespace {
  
  template <class T> struct param {

    //The integrand function
    const T & function;
    
  };

  template<class T> double integrand(double x , void * p) {
    //Units of the argument and return type
    const typename T::ValType ValUnit = 
      TypeTraits<typename T::ValType>::baseunit();
    const typename T::ArgType ArgUnit = 
      TypeTraits<typename T::ArgType>::baseunit();    

    const T & f = ((struct param<T> *)p)->function;
    return f(x * ArgUnit ) / ValUnit;
  }
  
}

namespace Herwig {
using namespace ThePEG;


template <class T>
inline typename BinaryOpTraits<typename T::ValType,
			       typename T::ArgType>::MulT
GSLIntegrator::value(const T & fn, 
		     const typename T::ArgType lower, 
		     const typename T::ArgType upper) const {
  typename BinaryOpTraits<typename T::ValType,
			       typename T::ArgType>::MulT error;
  return value(fn,lower,upper,error);
}


template <class T>
inline typename BinaryOpTraits<typename T::ValType,
			       typename T::ArgType>::MulT
GSLIntegrator::value(const T & fn, 
		     const typename T::ArgType lower, 
		     const typename T::ArgType upper,
		     typename BinaryOpTraits<typename T::ValType,
		     typename T::ArgType>::MulT & error) const {
  typedef typename T::ValType ValType;
  typedef typename T::ArgType ArgType;
  const ValType ValUnit = TypeTraits<ValType>::baseunit();
  const ArgType ArgUnit = TypeTraits<ArgType>::baseunit();
  
  double result(0.), error2(0.);
  
  param<T> parameters = { fn };
  gsl_function integrationFunction;
  integrationFunction.function = &integrand<T>;
  integrationFunction.params = &parameters;

  gsl_integration_workspace * workspace = 
    gsl_integration_workspace_alloc(_nbins);
  //do integration
  //Want to check error messages ourselves
  gsl_error_handler_t * oldhandler = gsl_set_error_handler_off();
  int status = gsl_integration_qags(&integrationFunction, lower/ArgUnit, 
				    upper/ArgUnit, _abserr, _relerr, _nbins, 
				    workspace, &result, &error2);
  if( status > 0 ) {
    CurrentGenerator::log() << "An error occurred in the GSL "
      "integration subroutine:\n";
    switch( status ) {
    case GSL_EMAXITER: 
      CurrentGenerator::log() << "The maximum number of subdivisions "
	"was exceeded.\n";
      break;
    case GSL_EROUND: 
      CurrentGenerator::log() << "Cannot reach tolerance because of "
	"roundoff error, or roundoff error was detected in the "
	"extrapolation table.\n";
      break;
    case GSL_ESING:
      CurrentGenerator::log() << "A non-integrable singularity or "
	"other bad integrand behavior was found in the integration "
	"interval.\n";
      break;
    case GSL_EDIVERGE:
      CurrentGenerator::log() << "The integral is divergent, "
	"or too slowly convergent to be integrated numerically.\n"; 
      break;
    default:
      CurrentGenerator::log() << "A general error occurred with code " 
			      << status << '\n';
    }
    result = 0.;
  }
  gsl_set_error_handler(oldhandler);
  gsl_integration_workspace_free(workspace);

  //fix units and return
  error = error2* ValUnit * ArgUnit;
  return result * ValUnit * ArgUnit;
}

}
