// -*- C++ -*-
//
// YFSFormFactors.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the YFSFormFactors class.
//

#include "SusyLoopIntegral.h"
#include <cassert>

using namespace Herwig;

const double SusyLoopIntegral::epsi = 1e-8;

double SusyLoopIntegral::B0(Energy2 s, Energy m1, Energy m2, Energy2 mu2) {
  Energy2 m12(sqr(m1)),m22(sqr(m2));
  if(s==ZERO) {
    if(m12 == m22 ) {
      return -log(m12/mu2);
    }
    else {
      if ( m12 != ZERO && m22 != ZERO ) {
	return 1. - m12/(m12-m22)*log(m12/mu2)
	  + m22/(m12-m22)*log(m22/mu2);
      }
      else if( m12 == ZERO ) {
	return 1. - log(m22/mu2);
      }
      else if( m22 == ZERO ) {
	return 1. - log(m12/mu2);
      }
      else
	assert(false);
    }
  }
  else {
    if ( m12 == ZERO && m22 == ZERO ) {
      return 2. - log(s/mu2);
    }
    else if ( m12 == s  && m22 == ZERO ) {
      return 2. - log(m12/mu2);
    }
    else if ( m22 == s && m12 == ZERO ) { 
      return 2. - log(m22/mu2);
    }
    else if (m12 ==ZERO) {
      return 2. - (s-m22)/s*log( abs(m22-s)/m22 ) - log(m22/mu2);
    }
    else if (m22 == ZERO) {
      return 2. - (s-m12)/s*log( abs(m12-s)/m12 ) - log(m12/mu2);
    }
    else {
      complex<double> arg = double((sqr(s)+sqr(m12)+sqr(m22)
				    -2.*(s*m12+s*m22+m12*m22))*UnitRemoval::InvE4);
      complex<Energy2> zkappa = sqrt(arg)*UnitRemoval::E2;
      complex<double> x1 = 0.5*(s-m22+m12+zkappa)/s;
      complex<double> x2 = 0.5*(s-m22+m12-zkappa)/s;
      return real( 2. + log(mu2/m22) + x1*log(1.-1./x1) 
		   + x2*log(1.-1./x2));
    }
  }
}

InvEnergy2 SusyLoopIntegral::B0P(Energy2 s, Energy m1, Energy m2, Energy2 mu2) {
  Energy2 m12(sqr(m1)),m22(sqr(m2));
  if(s==ZERO) {
    if( m12 == m22 ) {
      return 1./6./m12;
    }
    else {
      return ( 0.5*(m12+m22) - m12*m22/(m12-m22)*log(m12/m22) )/sqr(m12-m22);
    }
  }
  else if( s == m12 && m22 == ZERO ) {
    return ( -1. + 0.5*log(m12/mu2) )/m12;
  }
  else if( s == m22 && m12 == ZERO ) {
    return ( -1. + 0.5*log(m22/mu2) )/m22;
  }
  else if( m12 == ZERO && m22 != ZERO ) {
    return ( -1. - m22/s*log(abs(m22-s)/m22) )/s;
  }
  else if( m22 == ZERO && m12 != ZERO ) {
    return ( -1. - m12/s*log(abs(m12-s)/m12) )/s;
  }
  else {
    complex<double> arg = double((sqr(s)+sqr(m12)+sqr(m22)
				  -2.*(s*m12+s*m22+m12*m22))*UnitRemoval::InvE4);
    complex<Energy2> zkappa = sqrt(arg)*UnitRemoval::E2;
    complex<double> x1 = 0.5*(s-m22+m12+zkappa)/s;
    complex<double> x2 = 0.5*(s-m22+m12-zkappa)/s;
    return real( -1. + ( x1*(1. - x1) * log(1. - 1./x1) -
			 x2*(1. - x2) * log(1. - 1./x2) )/(x1-x2) )/s;
  }
}

complex<InvEnergy2> SusyLoopIntegral::C0(Energy2 p1, Energy2 p2, Energy2 p3,
					 Energy m1, Energy m2, Energy m3) {
  static const Energy2 eps = 1e-8*GeV2;
  static const Complex im(0.,1.);
  static const complex<InvEnergy2> ieps(ZERO,1e-17/GeV2);  
  Energy2 q10 = p1 != ZERO ? p1 : eps;
  Energy2 q20 = p3 != ZERO ? p3 : eps;
  Energy2 q21 = p2 != ZERO ? p2 : eps;
  Energy2 sm0 = sqr(m1), sm1 = sqr(m2), sm2 = sqr(m3);
  Energy2 R[3] = {p2,p3,p1};
  Energy2 alpha = abs( Complex(kappa(q10,q21,q20)*UnitRemoval::InvE2 ))*UnitRemoval::E2;
  complex<Energy2 > alp[3]={kappa(q21,sm1,sm2)*(1.+Complex(ieps*q21)),
			    kappa(q20,sm2,sm0)*(1.+Complex(ieps*q20)),
			    kappa(q10,sm0,sm1)*(1.+Complex(ieps*q10))};
  Complex x[3][3];
  x[0][1] = 0.5*(q21 - sm1 + sm2 + alp[0])/q21;
  x[0][2] = 0.5*(q21 - sm1 + sm2 - alp[0])/q21;
  x[1][1] = 0.5*(q20 - sm2 + sm0 + alp[1])/q20;
  x[1][2] = 0.5*(q20 - sm2 + sm0 - alp[1])/q20;
  x[2][1] = 0.5*(q10 - sm0 + sm1 + alp[2])/q10;
  x[2][2] = 0.5*(q10 - sm0 + sm1 - alp[2])/q10;
  Complex y0[3];
  y0[0] = (q21*(q21-q20-q10+2*sm0-sm1-sm2) - (q20-q10)*(sm1-sm2)
	   + alpha*(q21-sm1+sm2))/2/alpha/q21;
  y0[1] = (q20*(q20-q10-q21+2*sm1-sm2-sm0) - (q10-q21)*(sm2-sm0)
	   + alpha*(q20-sm2+sm0))/2/alpha/q20;
  y0[2] = (q10*(q10-q21-q20+2*sm2-sm0-sm1) - (q21-q20)*(sm0-sm1)
	   + alpha*(q10-sm0+sm1))/2/alpha/q10;
  Complex y[3][3];
  y[0][1] = y0[0] - x[0][1];
  y[0][2] = y0[0] - x[0][2];
  y[1][1] = y0[1] - x[1][1];
  y[1][2] = y0[1] - x[1][2];
  y[2][1] = y0[2] - x[2][1];
  y[2][2] = y0[2] - x[2][2];
  Complex output(0.);
  for(unsigned int i=0;i<3;++i) {
    for(unsigned int j=1;j<3;++j) {
      output += spenceFunction((y0[i]-1.)/y[i][j]) - 
	spenceFunction(y0[i]/y[i][j]);
      Complex cx = eta(1.-x[i][j],1./y[i][j]);
      if(abs(cx) != 0.) output += cx*log((y0[i]-1.)/y[i][j]);
      Complex cy = eta(-x[i][j],1./y[i][j]);
      if(abs(cy) != 0.) output -= cy*log(y0[i]/y[i][j]);
    }
    Complex cx = eta(-x[i][1],-x[i][2]);
    if(abs(cx)!=0.) output -= cx*log((1.-y0[i])/(-y0[i]));
    Complex cy = eta(y[i][1],y[i][2]);
    if(abs(cy)!=0) output += cy*log((1.-y0[i])/(-y0[i]));
    if( R[i] < ZERO && (y[i][1]*y[i][2]).imag() < 0. ) {
      output += Constants::twopi*im*log((1.-y0[i])/(-y0[i]));
    }
  }
  return output/alpha;
}

complex<InvEnergy4> 
SusyLoopIntegral::D0(Energy2 p1,Energy2 p2,Energy2 p3,Energy2 p4,
		     Energy2 p12,Energy2 p23,
		     Energy m1,Energy m2,Energy m3,Energy m4){
  Energy mm1=m1,mm2=m2,mm3=m3,mm4=m4;
  Energy2 m12=sqr(m1),m22=sqr(m2),m32=sqr(m3),m42=sqr(m4);
  Energy2 q1=p1,q2=p2,q3=p3,q4=p4,q12=p12,q23=p23;
  // at least one zero mass
  if(mm1==ZERO || mm2==ZERO || mm3==ZERO || mm4==ZERO) {
    // permute until mm3 is zero
    while (true) {
      Energy mm0=mm1;
      mm1=mm2;
      mm2=mm3;
      mm3=mm4;
      mm4=mm0;
      Energy2 m02=m12;
      m12=m22;
      m22=m32;
      m32=m42;
      m42=m02;
      Energy2 q00=q12;
      q12=q23;
      q23=q00;
      Energy2 q0=q1;
      q1=q2;
      q2=q3;
      q3=q4;
      q4=q0;
      if(mm3 == ZERO && (mm1 !=ZERO || mm2 != ZERO ) && 
	 mm4 == ZERO) continue;
      if(mm3 == ZERO) break;
    };
    // only mm3 is zero
    if(mm1!=ZERO&&mm2!=ZERO &&mm4!=ZERO) {
      double eps(1e-17);
      Complex ieps(0.,eps);
      Energy m[5]={ZERO,mm1,mm2,10.*GeV,mm4};
      Energy2 p[5][5];
      p[1][2]=q1;
      p[1][3]=q12;
      p[1][4]=q4;
      p[2][3]=q2;
      p[2][4]=q23;
      p[3][4]=q3;
      double k[5][5];
      Complex r[5][5],rs[5][5];
      for(int j=2;j<5;++j) {
	for(int i=1;i<j;++i) {
	  k[i][j] = (sqr(m[i])+sqr(m[j])-p[i][j])/m[i]/m[j];
	  if(i==3) k[i][j] -= m[i]/m[j];
	  if(j==3) k[i][j] -= m[j]/m[i];
	  r[i][j] = quadraticSolution(1.,-k[i][j],1.);
	  if(r[i][j].imag()==0.) {
	    rs[i][j] = quadraticSolution(1.,Complex(-k[i][j],eps),1.);
	  }
	  else {
	    rs[i][j] = r[i][j];
	  }
	}
      }
      Complex ss[5]={0.,rs[1][2],rs[2][3],rs[3][4],rs[1][4]};
      Complex aa = k[3][4]/r[2][4]-k[2][3];
      Complex bb = k[1][3]*(1./r[2][4]-r[2][4])+k[1][2]*k[3][4]-k[1][4]*k[2][3];
      Complex cc = k[1][2]*k[1][3]-k[1][3]*k[1][4]*r[2][4]+r[2][4]*k[3][4]-k[2][3];
      Complex dd = k[2][3]-r[2][4]*k[3][4];
      Complex xx[3]={0.,quadraticSolution(aa,bb,cc+ieps*dd),0.};
      xx[2] = (cc+ieps*dd)/aa/xx[1];
      Complex x[3][5];
      for(int i=1;i<3;++i) {
	x[i][1] = xx[i]/r[2][4];
	x[i][2] = xx[i]/r[2][4]*r[1][3];
	x[i][3] = xx[i]*r[1][3];
	x[i][4] = xx[i];
      }
      Complex output(0.);
      for(int i=1;i<3;++i) {
	output += (2.*double(i)-3.)*
	  (spenceFunction(1.+ss[4]*x[i][4]) - 
	   spenceFunction(1.+ss[1]*x[i][1]) +
	   spenceFunction(1.+x[i][4]/ss[4]) -
	   spenceFunction(1.+x[i][1]/ss[1]) +
	   eta(-x[i][4],ss[4])*log(1.+ss[4]*x[i][4]) - 
	   eta(-x[i][1],ss[1])*log(1.+ss[1]*x[i][1]) +
	   eta(-x[i][4],1./ss[4])*log(1.+x[i][4]/ss[4]) - 
	   eta(-x[i][1],1./ss[1])*log(1.+x[i][1]/ss[1]) - 
	   spenceFunction(1.+x[i][4]*(k[3][4]-ieps)/(k[1][3]-ieps)) +
	   spenceFunction(1.+x[i][1]*(k[2][3]-ieps)/(k[1][3]-ieps)) - 
	   eta(-x[i][4],(k[3][4]-ieps)/(k[1][3]-ieps))*
	   log(1.+x[i][4]*(k[3][4]-ieps)/(k[1][3]-ieps)) +
	   eta(-x[i][1],(k[2][3]-ieps)/(k[1][3]-ieps))*
	   log(1.+x[i][1]*(k[2][3]-ieps)/(k[1][3]-ieps)));
	Complex h(0.);
	if(r[2][4].imag()!=0.) {
	  h = eta(-1./xx[i],r[2][4]);
	}
	else if(r[2][4].real()<0.) {
	  Complex hh = -1./xx[i];
	  double im1 = hh.imag();
	  double im2 = rs[2][4].imag();
	  if(im1>0.&&im2>0.) h = -Complex(0.,Constants::twopi);
	  if(im1<0.&&im2<0.) h =  Complex(0.,Constants::twopi);
	}
	output += (2.*double(i)-3.)*h*
	  ( log( (k[1][2]-r[2][4]*k[1][4] +
		  xx[i]*(1./r[2][4]-r[2][4]))/dd ) +
	    log(k[1][3]-ieps) );
      }
      return output/m[1]/m[2]/m[3]/m[4]/aa/(xx[1]-xx[2]);
    }
    // mm2 and mm3 are zero
    else if (mm1 != ZERO && mm4 !=ZERO && mm2 == ZERO) {
      double eps(1e-17);
      Complex ieps(0.,eps);
      Energy m[5]={ZERO,mm1,10.*GeV,10.*GeV,mm4};
      Energy2 p[5][5];
      p[1][2]=q1;
      p[1][3]=q12;
      p[1][4]=q4;
      p[2][3]=q2;
      p[2][4]=q23;
      p[3][4]=q3;
      double k[5][5];
      Complex r[5][5],rs[5][5];
      for(int j=2;j<5;++j) {
	for(int i=1;i<j;++i) {
	  k[i][j] = (sqr(m[i])+sqr(m[j])-p[i][j])/m[i]/m[j];
	  if(i==2) k[i][j] -= m[i]/m[j];
	  if(j==2) k[i][j] -= m[j]/m[i];
	  if(i==3) k[i][j] -= m[i]/m[j];
	  if(j==3) k[i][j] -= m[j]/m[i];
	  r[i][j] = quadraticSolution(1.,-k[i][j],1.);
	  if(r[i][j].imag()==0.) {
	    rs[i][j] = quadraticSolution(1.,Complex(-k[i][j],eps),1.);
	  }
	  else {
	    rs[i][j] = r[i][j];
	  }
	}
      }
      Complex ss[5]={0.,rs[1][2],rs[2][3],rs[3][4],rs[1][4]};
      Complex aa = k[2][4]*k[3][4] - k[2][3];
      Complex bb = k[1][3]*k[2][4] + k[1][2]*k[3][4]-k[1][4]*k[2][3];
      Complex cc = k[1][2]*k[1][3] - k[2][3];
      Complex dd = k[2][3];
      Complex xx[3] = {0.,quadraticSolution(aa,bb,cc+ieps*dd),0.};
      xx[2] = (cc+ieps*dd)/aa/xx[1];
      Complex x[3][5];
      for(int i=1;i<3;++i) {
	x[i][1] = xx[i]/r[2][4];
	x[i][2] = xx[i]/r[2][4]*r[1][3];
	x[i][3] = xx[i]*r[1][3];
	x[i][4] = xx[i];
      }
      Complex output(0.);
      for(int i=1;i<3;++i) {
	output += (2.*double(i)-3.)*
	  (spenceFunction(1.+ss[4]*x[i][4]) + 
	   spenceFunction(1.+x[i][4]/ss[4]) +
	   eta(-x[i][4],ss[4])*log(1.+ss[4]*x[i][4]) +
	   eta(-x[i][4],1./ss[4])*log(1.+x[i][4]/ss[4]) - 
	   spenceFunction(1.+xx[i]*(k[3][4]-ieps)/(k[1][3]-ieps)) - 
	   spenceFunction(1.+xx[i]*(k[2][4]-ieps)/(k[1][2]-ieps)) -
	   eta(-xx[i],(k[3][4]-ieps)/(k[1][3]-ieps))*
	   log(1.+xx[i]*(k[3][4]-ieps)/(k[1][3]-ieps)) - 
	   eta(-xx[i],(k[2][4]-ieps)/(k[1][2]-ieps))*
	   log(1.+xx[i]*(k[2][4]-ieps)/(k[1][2]-ieps)) +
	   log(-xx[i])*( log(k[1][2]-ieps) + log(k[1][3]-ieps)
			 - log(k[2][3]-ieps) ));
      }
      return output/m[1]/m[2]/m[3]/m[4]/aa/(xx[1]-xx[2]);
    }
    else
      assert(false);
  }
  // all masses zero zero
  double eps(1e-18);
  Complex ieps(0.,eps);
  Energy m[5];
  m[0] = ZERO;
  Energy2 p[5][5];
  if(abs((sqr(mm1)+sqr(mm3)-q12)/mm1/mm3)<2.) {
    m[1] = mm2;
    m[2] = mm3;
    m[3] = mm4;
    m[4] = mm1;
    p[1][2]=q2;
    p[1][3]=q23;
    p[1][4]=q1;
    p[2][3]=q3;
    p[2][4]=q12;
    p[3][4]=q4;
  }
  else {
    m[1] = mm1;
    m[2] = mm2;
    m[3] = mm3;
    m[4] = mm4;
    p[1][2]=q1;
    p[1][3]=q12;
    p[1][4]=q4;
    p[2][3]=q2;
    p[2][4]=q23;
    p[3][4]=q3;
  }
  double k[5][5];
  Complex r[5][5],rs[5][5];
  for(int j=2;j<5;++j) {
    for(int i=1;i<j;++i) {
      k[i][j] = (sqr(m[i])+sqr(m[j])-p[i][j])/m[i]/m[j];
      r[i][j] = quadraticSolution(1.,-k[i][j],1.);
      if(r[i][j].imag()==0.) {
	rs[i][j] = quadraticSolution(1.,Complex(-k[i][j],eps),1.);
      }
      else {
	rs[i][j] = r[i][j];
      }
    }
  }
  Complex ss[5]={0.,rs[1][2],rs[2][3],rs[3][4],rs[1][4]};
  Complex s0[5]={0.,r[1][2],r[2][3],r[3][4],r[1][4]};
  Complex aa = k[3][4]/r[2][4]+r[1][3]*k[1][2]-k[1][4]*r[1][3]/r[2][4]-k[2][3];
  Complex bb = (r[2][4]-1./r[2][4])*(r[1][3]-1./r[1][3])
    +k[1][2]*k[3][4]-k[1][4]*k[2][3];
  Complex cc = k[1][2]/r[1][3]+r[2][4]*k[3][4]-k[1][4]*r[2][4]/r[1][3]-k[2][3];
  Complex dd = k[2][3]-r[1][3]*k[1][2]-r[2][4]*k[3][4]+r[1][3]*r[2][4]*k[1][4];
  Complex xx[3]={0.,quadraticSolution(aa,bb,cc+ieps*dd),0.};
  xx[2] = (cc+ieps*dd)/aa/xx[1];
  Complex xx0[3]={0.,quadraticSolution(aa,bb,cc),0.};
  xx0[2] = cc/aa/xx0[1];
  if(abs(xx0[1]-xx[2]) < abs(xx0[1]-xx[1])) swap(xx0[1],xx0[2]);
  
  Complex g[3],x[3][5],x0[3][5];
  for(int i=1;i<3;++i) {
    g[i] = (aa*(xx[i]-xx[3-i])).real() >=0. ? 1. : -1.;
    x[i][1]  = xx[i]/r[2][4];
    x0[i][1] = xx0[i]/r[2][4];
    x[i][2]  = xx[i]/r[2][4]*r[1][3];
    x0[i][2] = xx0[i]/r[2][4]*r[1][3];
    x[i][3]  = xx[i]*r[1][3];
    x0[i][3] = xx0[i]*r[1][3];
    x[i][4]  = xx[i];
    x0[i][4] = xx0[i];
  }
  Complex output(0.);
  for(int i=1;i<3;++i) {
    for(int j=1;j<5;++j) {
      double sign = (x[i][j]/ss[j]).imag()>=0. ? 1. : -1.;
      Complex a1 = 1.+x0[i][j]*s0[j] + abs(1.+x0[i][j]*s0[j])*ieps*sign;
      Complex a2 = 1.+x0[i][j]/s0[j] + abs(1.+x0[i][j]/s0[j])*ieps*sign;
      output += pow(-1.,i+j)*
	(spenceFunction(a1)+eta(-x[i][j],ss[j])*log(a1) + 
	 spenceFunction(a2)+eta(-x[i][j],1./ss[j])*log(a2));
    }
  }
  if(r[1][3].imag() == 0. ) {
    for(int i=1;i<3;++i) {
      Complex a1 = (k[1][3]-2.*r[1][3])/xx0[i]-r[1][3]*k[1][4]+k[3][4];
      Complex a2 = ((k[2][4]-2.*r[2][4])*r[1][3]*xx0[i]
		    -r[2][4]*k[3][4]+k[2][3])/dd;
      Complex a3 = (k[1][3]-2.*r[1][3])*r[2][4]/xx0[i]-r[1][3]*k[1][2]+k[2][3];
      Complex a4 = ((k[2][4]-2.*r[2][4])*xx0[i]-r[2][4]*k[1][4]+k[1][2])/dd;
      Complex l1 = log( a1-abs(a1)*ieps );
      double sign = r[1][3].real()*rs[2][4].imag() >=0. ? 1. : -1.;
      Complex l2 = log( a2+abs(a2)*ieps*g[i]*sign );
      Complex l3 = log( a3-abs(a3)*ieps );
      sign = rs[2][4].imag() >= 0. ? 1. : -1.;
      Complex l4 = log( a4+abs(a4)*ieps*g[i]*sign);
      output += (3.-2.*double(i))*
	( etas( -xx[i],r[1][3],rs[1][3] )*
	  ( log(r[1][3]*xx[i]) + l1 + l2 ) +
	  etas( -xx[i],1./r[2][4],1./rs[2][4] )*
	  ( log(xx[i]/r[2][4]) + l3 + l4 ) - 
	  ( etas( -xx[i],r[1][3]/r[2][4],rs[1][3]/rs[2][4] ) +
	    eta ( rs[1][3],1./rs[2][4] ))*
	  ( log(xx[i]*r[1][3]/r[2][4]) + l3 + l2 ) +
	  eta( rs[1][3],1./rs[2][4] )*
	  etas(-xx[i],-r[1][3]/r[2][4],-rs[1][3]/rs[2][4])   
	  );
    }
  }
  else {
    for(int i=1;i<3;++i) {
      Complex l1 = log( r[2][4]/xx0[i]+xx0[i]/r[2][4]+k[1][2]
			-xx0[i]/r[2][4]*eps*bb*g[i] );
      Complex l2 = log( r[1][3]*xx0[i]+1./xx0[i]/r[1][3]+k[3][4]
			-xx0[i]*r[1][3]*eps*bb*g[i] );
      Complex l3 = log( r[1][3]/r[2][4]*xx0[i]+r[2][4]/xx0[i]/r[1][3]+k[2][3]
			-xx0[i]*r[1][3]/r[2][4]*eps*bb*g[i] );
      double sign = bb.real() >0. ? 1. : -1.;
      output += (3.-2.*double(i))*
	( eta(-xx[i],1./r[2][4])*( log(xx[i]/r[2][4]) + l1 ) + 
	  eta(-xx[i],   r[1][3])*( log(r[1][3]*xx[i]) + l2 ) -
	  ( eta(-xx[i],r[1][3]/r[2][4])+eta(r[1][3],1./r[2][4]) )*
	  ( log(xx[i]*r[1][3]/r[2][4]) + l3 ) +
	  eta(r[1][3],1./r[2][4])*eta(-xx[i],-r[1][3]/r[2][4])*
	  (1.-g[i]*sign));
    }
  }
  return output/m[1]/m[2]/m[3]/m[4]/aa/(xx[1]-xx[2]);
}

Complex SusyLoopIntegral::spenceFunction(Complex z) {
  double b[10] = {0.,
		  0.1666666666666666666666666667,
		  -0.0333333333333333333333333333,
		  0.0238095238095238095238095238,
		  -0.0333333333333333333333333333,
		  0.0757575757575757575757575758,
		  -0.2531135531135531135531135531,
		  1.1666666666666666666666666667,
		  -7.09215686274509804         ,
		  54.97117794486215539};
  double rz = z.real();
  double az = abs(z);
  double a1 = abs(1.-z);
  if(az<1e-20) return -log(1.-z);
  if( abs(rz-1.) < 1e-18 && abs(z.imag()) < 1e-18)
    return 1.64493406684822643;
  if(rz > 0.5 ) {
    if(a1>1.) {
      Complex w = log(1.-1./z);
      Complex sum = w - 0.25*sqr(w);
      Complex u = w;
      if(abs(u)>=1e-10) {
	for(int k=1;k<10;++k) {
	  u *= sqr(w)/Complex(2*k*(2*k+1));
	  if(abs(u*b[k]/sum)<1e-20) break;
	  sum += u*b[k];
	}     
      }
      return sum+3.28986813369645287+.5*sqr(log(z-1.))-log(z)*log(1.-z);
    }
    else {
      Complex w = -log(z);
      Complex sum = w - 0.25*sqr(w);
      Complex u = w;
      if(abs(u)>=1e-10) {
	for(int k=1;k<10;++k) {
	  u *= sqr(w)/Complex(2*k*(2*k+1));
	  if(abs(u*b[k]/sum)<1e-20) break;
	  sum += u*b[k];
	}
      }
      return -sum+1.64493406684822643-log(z)*log(1.-z);
    }     
  }
  else if(az > 1.) {
    Complex w = -log(1.-1./z);
    Complex sum = w-0.25*sqr(w);
    Complex u = w;
    if(abs(u)>=1e-10) {                               
      for(int k=1;k<10;++k) {
	u *= sqr(w)/Complex(2*k*(2*k+1));
	if(abs(u*b[k]/sum)<1e-20) break;
	sum += u*b[k];
      }
    }
    return -sum-1.64493406684822643-.5*sqr(log(-z));
  }
  else {
    Complex w = -log(1.-z);
    Complex sum = w - 0.25*sqr(w);
    Complex u = w;
    if(abs(u)>=1e-10) {                               
      for(int k=1;k<10;++k) {
	u *= sqr(w)/Complex(2*k*(2*k+1));
	if(abs(u*b[k]/sum)<1e-20) break;
	sum += u*b[k];
      }
    }
    return sum;
  }
}                
