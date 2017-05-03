// -*- C++ -*-
//
// Maths.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include <iostream>
#include "Maths.h"
namespace Herwig {
namespace Math {
using ThePEG::Complex;
using std::cout;
using std::endl;

namespace {
  inline Complex Li2Prod(Complex y,Complex y2)
  {
    static const double a1   =-0.250000000000000e0,a2   =-0.111111111111111e0;
    static const double a3   =-0.010000000000000e0,a4   =-0.017006802721088e0;
    static const double a5   =-0.019444444444444e0,a6   =-0.020661157024793e0;
    static const double a7   =-0.021417300648069e0,a8   =-0.021948866377231e0;
    static const double a9   =-0.022349233811171e0,a10  =-0.022663689135191e0;
    return y*(1.+a1*y*(1.+a2*y*(1.+a3*y2*(1.+a4*y2*(1.+a5*y2*
	   (1.+a6*y2*(1.+a7*y2*(1.+a8*y2*(1.+a9*y2*(1.+a10*y2))))))))));
  }

inline long double Li2Prod(long double y,long double y2)
{
  static const long double a1   =-0.250000000000000e0,a2   =-0.111111111111111e0;
  static const long double a3   =-0.010000000000000e0,a4   =-0.017006802721088e0;
  static const long double a5   =-0.019444444444444e0,a6   =-0.020661157024793e0;
  static const long double a7   =-0.021417300648069e0,a8   =-0.021948866377231e0;
  static const long double a9   =-0.022349233811171e0,a10  =-0.022663689135191e0;
  return y*(1.+a1*y*(1.+a2*y*(1.+a3*y2*(1.+a4*y2*(1.+a5*y2*
	   (1.+a6*y2*(1.+a7*y2*(1.+a8*y2*(1.+a9*y2*(1.+a10*y2))))))))));
  }
}


Complex Li2(Complex x)
{
  Complex z;
  static double zeta2= 1.644934066848226e0;
  double xr(real(x)),xi(imag(x)),r2(xr*xr+xi*xi);
  if(r2>1.&&xr/r2>0.5)
    {
      z=-log(1./x);
      return Li2Prod(z,z*z)+zeta2-log(x)*log(1.-x)+0.5*log(x)*log(x);
    }
  else if (r2>1.&&(xr/r2)<=0.5) 
    {
      z=-log(1.-1./x);
    return -Li2Prod(z,z*z)-zeta2-0.5*pow(log(-x),2);
    }
  else if(r2==1.&&xi==0.)
    {
      if(xr>0){return zeta2;}
      else{return -0.5*zeta2;}
    }
  else if(r2<=1.&&xr>0.5)
    {
      z=-log(x);
      return -Li2Prod(z,z*z)+zeta2-log(x)*log(1.-x);
    }
  else
    {
      z=-log(1.-x);
      return Li2Prod(z,z*z);
    }
}

long double ReLi2(long double x)
{
  long double z;
  static long double zeta2= 1.644934066848226e0;
  long double output;
  if(x>1.&&x<2.)
    {
      z=-log(1./x);
      output= Li2Prod(z,z*z)+zeta2-log(x)*log(x-1.)+0.5*log(x)*log(x);
    }
  else if (x>1.||x<-1.) 
    {
      z=-log(1.-1./x);
      if(x<-1){output= -Li2Prod(z,z*z)-zeta2-0.5*pow(log(-x),2);}
      else{output= -Li2Prod(z,z*z)-0.5*pow(log(x),2)+2.*zeta2;}
    }
  else if(x== 1.){output= zeta2;}
  else if(x==-1.){output= -0.5*zeta2;}
  else if(x>0.5)
    {
      z=-log(x);
      output= -Li2Prod(z,z*z)+zeta2-log(x)*log(1.-x);
    }
  else
    {
      z=-log(1.-x);
      output= Li2Prod(z,z*z);
    }
  return output;
}



}
}
