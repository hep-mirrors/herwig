// -*- C++ -*-

#include "Math.h"
namespace Herwig {
namespace Math {
using ThePEG::Complex;

Complex Li2(Complex x)
{
  Complex z;
  static double zeta2= 1.644934066848226e0;
  double xr(real(x)),xi(imag(x)),r2(xr*xr+xi*xi);
  if(r2>1.&&xr/r2>0.5)
    {
      z=-log(1./x);
      return Li2Prod(z,z*z)+zeta2-log(x)*log(1.-x)+0.5*pow(log(x),2);
    }
  else if (r2>1.&&(xr/r2)<=0.5) 
    {
      z=-log(1.-1./x);
    return -Li2Prod(z,z*z)-zeta2-0.5*pow(log(-x),2);
    }
  else if(r2==1.&&xi==0.)
    {
      return zeta2;
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

}
}
