// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GaussianIntegral class.
//
//  Author: Peter Richardson
//

#include "GaussianIntegral.h"
#include <cmath>
#include <iostream>

namespace Herwig 
{
  
using namespace Genfun;
using std::endl;
  
GaussianIntegral::~GaussianIntegral() {} 

// the member that does the actual calculation
double GaussianIntegral::operator [] (const AbsFunction & function) const 
{
  // vector for the limits of the bin
  vector<double> lowerlim,upperlim;
  // start with the whole interval as 1 bin
  lowerlim.push_back(_lower);upperlim.push_back(_upper);
  // set the minimum bin width
  double xmin=_binwidth*std::abs(_upper-_lower);
  // counters for the number of function evals
  int neval=0;
  // and number of bad intervals
  int nbad=0;
  // the output value
  double output=0.;
  // the loop for the evaluation
  double mid,wid; unsigned int ibin,ix=0,iorder;
  double testvalue,value,tolerance;
  do 
    {
      // the bin we are doing (always the last one in the list)
      ibin = lowerlim.size()-1;
      // midpoint and width of the bin
      mid=0.5*(upperlim[ibin]+lowerlim[ibin]);
      wid=0.5*(upperlim[ibin]-lowerlim[ibin]);
      value=0.;
      iorder=0;
      // compute a trail value using sixth order GQ
      for(ix=0;ix<_weights[0].size();++ix)
	{
	  value+=_weights[0][ix]*(+function(mid+wid*_abscissae[0][ix])
				  +function(mid-wid*_abscissae[0][ix]));
	  ++neval;
	  if(neval>_maxeval)
	    {std::cerr << "Error in Gaussian Integral: Setting to zero" 
		       << std::endl;}
	}
      value *=wid;
      // compute more accurate answers using higher order GQ
      do 
	{
	  // use the next order of quadrature
	  testvalue=value;
	  ++iorder;
	  value=0.;
	  for(ix=0;ix<_weights[iorder].size();++ix)
	    {
	      value+=_weights[iorder][ix]*
		(+function(mid+wid*_abscissae[iorder][ix])
		 +function(mid-wid*_abscissae[iorder][ix]));
	      ++neval;
	      if(neval>_maxeval)
		{std::cerr << "Error in Gaussian Integral: Setting to zero" 
			   << std::endl;}
	    }
	  value *=wid;
	  tolerance=std::max(_abserr,_relerr*std::abs(value));
	}
      // keep going if possible and not accurate enough
      while(iorder<_weights.size()-1&&std::abs(testvalue-value)>tolerance);
      // now decide what to do
      // accept this value
      if(std::abs(testvalue-value)<tolerance)
	{
	  output+=value;
	  lowerlim.pop_back();upperlim.pop_back();
	}
      // bin too small to redivide contribution set to zero
      else if(wid<xmin)
	{
	  ++nbad;
	  lowerlim.pop_back(); upperlim.pop_back();
	}
      // otherwise split the bin into two
      else
	{
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
    {std::cerr << "Error in GaussianIntegral: Bad Convergence for " 
	       << nbad << "intervals" << std::endl;}
  // return the answer
  return output;
}

}
