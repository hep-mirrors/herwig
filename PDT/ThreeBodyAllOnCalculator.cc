// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyAllOnCalculator class.
//

#include "ThreeBodyAllOnCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeBodyAllOnCalculator.tcc"
#endif

namespace Herwig {
using namespace ThePEG;
using namespace Genfun;

ThreeBodyAllOnCalculator::~ThreeBodyAllOnCalculator() 
{
  delete _theInnerIntegrand;
  delete _theOuterIntegrand;
  //delete _theME;
  delete _Integrator;
}

NoPIOClassDescription<ThreeBodyAllOnCalculator> ThreeBodyAllOnCalculator::initThreeBodyAllOnCalculator;
// Definition of the static class description member.

void ThreeBodyAllOnCalculator::Init() {

  static ClassDocumentation<ThreeBodyAllOnCalculator> documentation
    ("The \\classname{ThreeBodyAllOnCalculator} class integrators a 1-3 matrix"
     " element to give the width");

}

// shift the variables for the outer integrand and give limits for the inner one
void ThreeBodyAllOnCalculator::outerVariables(const double & x, double & low, double & upp)
{
  // first convert the value of x into the value of souter
  if(_channelmass[_thechannel]>0)
    {_souter = _channelmass[_thechannel]*(_channelmass[_thechannel]+
					   _channelwidth[_thechannel]*tan(x));}
  else
    {
      _souter = pow(x,1./(_channelwidth[_thechannel]+1.));
    }
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
double ThreeBodyAllOnCalculator::innerIntegrand(const double & y)
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
      if(_channelmass[ix]>0)
	{
	  rm2=_channelmass[ix]*_channelmass[ix];
	  rw2 = _channelwidth[ix]*_channelwidth[ix];
	  term = (sjac-rm2)*(sjac-rm2)+rw2*rm2;
	  term = _channelweights[ix]*_channelmass[ix]*_channelwidth[ix]/term;
	}
      else
	{
	  term = _channelweights[ix]*(_channelwidth[ix]+1.)*
	    pow(sjac,_channelwidth[ix]);
	}
      jacdem+=term;
    }
  // now computer the matrix element
  Genfun::Argument a(7); a[0]=_m2[0];a[1]=s12;a[2]=s13;a[3]=s23;
  a[4]=_m[1];a[5]=_m[2];a[6]=_m[3];
  return (*_theME)(a)/jacdem;
}
// calculate the width for a given mass
Energy ThreeBodyAllOnCalculator::partialWidth(Energy2 q2) const
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
      if(_channelmass[ix]>0)
	{
	  upp = atan((upp-_channelmass[ix]*_channelmass[ix])/
		     _channelmass[ix]/_channelwidth[ix]);
	  low =  atan((low-_channelmass[ix]*_channelmass[ix])/
		      _channelmass[ix]/_channelwidth[ix]);
	}
      else
	{
	  upp = pow(upp,_channelwidth[ix]+1.);
	  low = pow(low,_channelwidth[ix]+1.);
	}
      // perform the integral using my Gaussian quadature class
      _thechannel=ix;
      GaussianIntegral *intb= new GaussianIntegral(low,upp);
      value = _channelweights[ix]*(*intb)[*_theOuterIntegrand];
      delete intb;
      sum+= value;
    }
  // final factors
  double fact=2.*pi*_m[0];
  fact*=fact*fact;
  sum =sum/fact/32.;
  return sum;
}
}

namespace Herwig {
using namespace Genfun;

// outer integral
FUNCTION_OBJECT_IMP(ThreeBodyAllOnOuter)
  
ThreeBodyAllOnOuter::ThreeBodyAllOnOuter(ThreeBodyAllOnCalculatorPtr in)
{
  _integrand=in;
  _Integrator= new GaussianIntegral(0.,1.);
}
  
ThreeBodyAllOnOuter::~ThreeBodyAllOnOuter() 
{delete _Integrator;}
  
ThreeBodyAllOnOuter::ThreeBodyAllOnOuter(const ThreeBodyAllOnOuter & right) 
  {  }
  

double ThreeBodyAllOnOuter::operator() (double x) const 
{
  double low,upp;
  _integrand->outerVariables(x,low,upp);
  _Integrator->resetLimits(low,upp);
  return (*_Integrator)[*(_integrand->InnerIntegrand())];
}

// inner integral
FUNCTION_OBJECT_IMP(ThreeBodyAllOnInner)
    
ThreeBodyAllOnInner::ThreeBodyAllOnInner(ThreeBodyAllOnCalculatorPtr in)
{_integrand=in;}

ThreeBodyAllOnInner::~ThreeBodyAllOnInner() {}
  
ThreeBodyAllOnInner::ThreeBodyAllOnInner(const ThreeBodyAllOnInner & right) 
{  }
  
double ThreeBodyAllOnInner::operator() (double x) const 
{
  // calculate the integrand
  double output=_integrand->innerIntegrand(x);
  return output;
}

}
