// -*- C++ -*-
//
// GaussianIntegrator.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the GaussianIntegrator class.
//

namespace Herwig {
using namespace ThePEG;

template <class T>
inline GaussianIntegrator::ValT<T>
GaussianIntegrator::value(const T & function, 
			  const typename T::ArgType lower, 
			  const typename T::ArgType upper) const {
  typedef typename T::ValType ValType;
  typedef typename T::ArgType ArgType;
  const ValType ValUnit = TypeTraits<ValType>::baseunit();
  const ArgType ArgUnit = TypeTraits<ArgType>::baseunit();

  // vector for the limits of the bin
  vector<double> lowerlim,upperlim;
  // start with the whole interval as 1 bin
  lowerlim.push_back(lower/ArgUnit);upperlim.push_back(upper/ArgUnit);
  // set the minimum bin width
  double xmin=_binwidth*abs(upper-lower)/ArgUnit;
  // counters for the number of function evals
  int neval=0;
  // and number of bad intervals
  int nbad=0;
  // the output value
  double output=0.;
  // the loop for the evaluation
  double mid,wid; unsigned int ibin,ix=0,iorder;
  double testvalue,value,tolerance;
  do {
    // the bin we are doing (always the last one in the list)
    ibin = lowerlim.size()-1;
    // midpoint and width of the bin
    mid=0.5*(upperlim[ibin]+lowerlim[ibin]);
    wid=0.5*(upperlim[ibin]-lowerlim[ibin]);
    value=0.;
    iorder=0;
    // compute a trail value using sixth order GQ
    for(ix=0;ix<_weights[0].size();++ix) {
      value+=_weights[0][ix]
	*( function((mid+wid*_abscissae[0][ix])*ArgUnit)
	  +function((mid-wid*_abscissae[0][ix])*ArgUnit)
	   )/ValUnit;
      ++neval;
      if(neval>_maxeval) 
	CurrentGenerator::log() << "Error in Gaussian Integrator: Setting to zero" 
				<< endl;
    }
    value *=wid;
    // compute more accurate answers using higher order GQ
    do {
      // use the next order of quadrature
      testvalue=value;
      ++iorder;
      value=0.;
      for(ix=0;ix<_weights[iorder].size();++ix) {
	value+=_weights[iorder][ix]*
	  ( function((mid+wid*_abscissae[iorder][ix])*ArgUnit)
	   +function((mid-wid*_abscissae[iorder][ix])*ArgUnit)
	    )/ValUnit;
	++neval;
	if(neval>_maxeval)
	   CurrentGenerator::log() << "Error in Gaussian Integrator: Setting to zero" 
				   << endl;
      }
      value *=wid;
      tolerance=max(_abserr,_relerr*abs(value));
    }
    // keep going if possible and not accurate enough
    while(iorder<_weights.size()-1&&abs(testvalue-value)>tolerance);
    // now decide what to do
    // accept this value
    if(abs(testvalue-value)<tolerance) {
      output+=value;
      lowerlim.pop_back();upperlim.pop_back();
    }
    // bin too small to redivide contribution set to zero
    else if(wid<xmin) {
      ++nbad;
      lowerlim.pop_back(); upperlim.pop_back();
    }
    // otherwise split the bin into two
    else {
      // reset the limits for the bin
      upperlim[ibin]=mid;
      // set up a new bin
      lowerlim.push_back(mid);
      upperlim.push_back(mid+wid);
    }
  }
  // keep going if there's still some bins to evaluate
  while(lowerlim.size()>0);
  // output an error message if needed
  if(nbad!=0)
    CurrentGenerator::log() << "Error in GaussianIntegrator: Bad Convergence for " 
			    << nbad << "intervals" << endl;
  // return the answer
  return output * ValUnit * ArgUnit;
}

}
