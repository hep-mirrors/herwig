// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyIntegrator class.
//
//  Author: Peter Richardson
//

#include "ThreeBodyIntegrator.h"
#include "ThePEG/Interface/ClassDocumentation.h"



namespace Herwig {

using namespace ThePEG;
using namespace Genfun;
  
ThreeBodyIntegrator::~ThreeBodyIntegrator() {}
  
void ThreeBodyIntegrator::Init() {
  
static ClassDocumentation<ThreeBodyIntegrator> documentation
  ("There is no documentation for the \\classname{ThreeBodyIntegrator} class");
 
}

// shift the variables for the outer integrand and give limits for the inner one
void ThreeBodyIntegrator::outerVariables(const double & x, double & low, double & upp)
{
  // first convert the value of x into the value of souter
  _souter = _channelmass[_thechannel]*(_channelmass[_thechannel]+
				       _channelwidth[_thechannel]*tan(x));
  // now the limits of the inner integral
  Energy ea(0.),eb(0.),eam(0.),ebm(0.);
  Energy rs=sqrt(_souter);
  switch(_channeltype[_thechannel])
    {
    case 1:
      ea = 0.5*(_souter-_m2[1]+_m2[2])/rs; eam=sqrt(ea*ea-_m2[2]);
      eb = 0.5*(_m2[0]-_souter-_m2[3])/rs; ebm=sqrt(eb*eb-_m2[3]);
      break;
    case 2:
      ea = 0.5*(_souter-_m2[1]+_m2[3])/rs; eam=sqrt(ea*ea-_m2[3]);
      eb = 0.5*(_m2[0]-_souter-_m2[2])/rs; ebm=sqrt(eb*eb-_m2[2]);
      break;
    case 3:
      ea = 0.5*(_souter-_m2[2]+_m2[3])/rs; eam=sqrt(ea*ea-_m2[3]);
      eb = 0.5*(_m2[0]-_souter-_m2[1])/rs; ebm=sqrt(eb*eb-_m2[1]);
      break;
    }
  Energy2 sum=(ea+eb)*(ea+eb);
  // calculate the limits
  low = sum-(eam+ebm)*(eam+ebm);
  upp = sum-(eam-ebm)*(eam-ebm);
}

// the integrand for the inner integral
double ThreeBodyIntegrator::innerIntegrand(const double & y)
{
  // set up the values of the s variables
  Energy2 s12(0.),s23(0.),s13(0.),m2sum=_m2[0]+_m2[1]+_m2[2]+_m2[3];
  switch(_channeltype[_thechannel])
    {
    case 1:
      s12 = _souter;
      s23 = y;
      s13 = m2sum-s12-s23;
      break;
    case 2:
      s23 = y;
      s13 = _souter;
      s12 = m2sum-s23-s13;
      break;
    case 3:
      s23 = _souter;
      s13 = y;
      s12 = m2sum-s23-s13;
      break;
    }
  // compute the jacobian
  // computer the denominator for the jacobian
  double jacdem=0.,term; Energy sjac(0.); Energy rm2,rw2;
  for(unsigned int ix=0,N=_channeltype.size();ix<N;++ix)
    {
      switch(_channeltype[ix])
	{
	case 1:
	  sjac = s12;
	  break;
	case 2:
	  sjac=s13;
	  break;
	case 3:
	  sjac=s23;
	  break;
	}
      rm2=_channelmass[ix]*_channelmass[ix];
      rw2 = _channelwidth[ix]*_channelwidth[ix];
      term = (sjac-rm2)*(sjac-rm2)+rw2*rm2;
      term = _channelweights[ix]/term;
      jacdem+=term;
    }
  // now computer the matrix element
  Genfun::Argument a(4); a[0]=_m2[0];a[1]=s12;a[2]=s13;a[3]=s23;
  return (*_theME)(a)/jacdem;
}

// calculate the width for a given mass
Energy ThreeBodyIntegrator::width(Energy2 q2)
{
  _m[0] = sqrt(q2);
  _m2[0]=q2;
  // check the decay is kinematically allowed
  if(_m[0]<_m[1]+_m[2]+_m[3]){return 0.;}
  // perform the integrals for all the different channels
  double upp(0.),low(0.),sum(0.),value;
  for(unsigned int ix=0,N=_channeltype.size();ix<N;++ix)
    {
      // work out the kinematic limits
      switch(_channeltype[ix])
	{
	case 1:
	  upp = (_m[0]-_m[3])*(_m[0]-_m[3]);
	  low = (_m[1]+_m[2])*(_m[1]+_m[2]);
	  break;
	case 2:
	  upp = (_m[0]-_m[2])*(_m[0]-_m[2]);
	  low = (_m[1]+_m[3])*(_m[1]+_m[3]);
	  break;
	case 3:
	  upp = (_m[0]-_m[1])*(_m[0]-_m[1]);
	  low = (_m[2]+_m[3])*(_m[2]+_m[3]);
	  break;
	}
      // transform them
      upp = atan((upp-_channelmass[ix]*_channelmass[ix])/
		 _channelmass[ix]/_channelwidth[ix]);
      low =  atan((low-_channelmass[ix]*_channelmass[ix])/
		  _channelmass[ix]/_channelwidth[ix]);
      // perform the integral using my Gaussian quadature class
      _thechannel=ix;
      GaussianIntegral *intb= new GaussianIntegral(low,upp);
      value = _channelweights[ix]*(*intb)[*_theOuterIntegrand];
      value=value/_channelmass[ix]/_channelwidth[ix];
      sum+= value;
    }
  // final factors
  double fact=2.*acos(-1.0)*_m[0];
  fact*=fact*fact;
  sum =sum/fact/32.;
  return sum;
}
}

namespace Herwig {

using namespace Genfun;

// outer integral
FUNCTION_OBJECT_IMP(ThreeBodyOuterIntegrand)
  
ThreeBodyOuterIntegrand::ThreeBodyOuterIntegrand(ThreeBodyIntegrator * in)
{_theIntegrator=in;}
  
ThreeBodyOuterIntegrand::~ThreeBodyOuterIntegrand() {}
  
ThreeBodyOuterIntegrand::ThreeBodyOuterIntegrand(const ThreeBodyOuterIntegrand & right) 
  {  }
  

double ThreeBodyOuterIntegrand::operator() (double x) const 
{
  // calculate the limits of the inner integral
  double low,upp;
  _theIntegrator->outerVariables(x,low,upp);
  // setup the inner integral
  GaussianIntegral *inte=new GaussianIntegral(low,upp);
  double output=(*inte)[*(_theIntegrator->InnerIntegrand())];
  delete inte;
  return output;
}

// inner integral
FUNCTION_OBJECT_IMP(ThreeBodyInnerIntegrand)
    
ThreeBodyInnerIntegrand::ThreeBodyInnerIntegrand(ThreeBodyIntegrator * in)
{_theIntegrator=in;}

ThreeBodyInnerIntegrand::~ThreeBodyInnerIntegrand() {}
  
ThreeBodyInnerIntegrand::ThreeBodyInnerIntegrand(const ThreeBodyInnerIntegrand & right) 
{  }
  
double ThreeBodyInnerIntegrand::operator() (double x) const 
{
  // calculate the integrand
  return _theIntegrator->innerIntegrand(x);
}

}
