// -*- C++ -*-
//
// GSLIntegrator.tcc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
      TypeTraits<typename T::ValType>::baseunit;
    const typename T::ArgType ArgUnit = 
      TypeTraits<typename T::ArgType>::baseunit;    

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
  typedef typename T::ValType ValType;
  typedef typename T::ArgType ArgType;
  const ValType ValUnit = TypeTraits<ValType>::baseunit;
  const ArgType ArgUnit = TypeTraits<ArgType>::baseunit;
  
  double result(0.), error(0.);
  
  param<T> parameters = { fn };
  gsl_function integrationFunction;
  integrationFunction.function = &integrand<T>;
  integrationFunction.params = &parameters;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(_nbins);
  //do integration
  gsl_integration_qags(&integrationFunction, lower/ArgUnit, upper/ArgUnit, 
		       _abserr, _relerr, _nbins, w, &result, &error); 
  gsl_integration_workspace_free(w);

  //fix units and return
  return result * ValUnit * ArgUnit;
}

}
