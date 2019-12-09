// -*- C++ -*-
//
// YFSFormFactors.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the YFSFormFactors class.
//

#include "YFSFormFactors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <cassert>

using namespace Herwig;

using Constants::pi;
using Herwig::Math::ReLi2;

const double  YFSFormFactors::_alpha=1./137.03599911;

const Energy  YFSFormFactors::_mgamma=1e-10*MeV;

const Energy2 YFSFormFactors::_tcut=1.e-11*GeV2;

const Energy YFSFormFactors::_ecut=1e-6*GeV;

double YFSFormFactors::ReBIF(Energy  m0      ,Energy  m1      , Energy2 t       ,
			     double  charge  ,bool    includegamma,
			     Energy  mgamma) {
  // mass squared for speed
  Energy2 m02(m0*m0),m12(m1*m1),nu(0.5*(m02+m12-t)),mprod(m0*m1);
  double Anu,vfinite;
  double output;
  // t>0
  if(t>_tcut) {
    // parameters
    Energy2 lambda(sqrt((nu-mprod)*(nu+mprod)));
    double eta(0.5*m12*t/lambda/(lambda+nu-m12)),zeta((lambda+nu)*eta/m12);
    // simple A functions for virtual piece
    InvEnergy2 A;
    if(lambda>1e-6*GeV2){A=(log((lambda+nu)/mprod)/lambda);}
    else{A=1./mprod;}
    double A1((m02-m12)/t*log(m0/m1)-2.*sqr(lambda)/t*A-2.);
    InvEnergy2 A3(A*log(2.*lambda/mprod)
		  +1./lambda*
		  (+0.25*(log((lambda+nu)/m02)+2.*log((lambda-nu+m02)/t  ))*
		   log((lambda+nu)/m02)
		   +0.25*(log((lambda+nu)/m12)-2.*log((lambda+nu-m12)/m12))*
		   log((lambda+nu)/m12)
		   +0.5*(log(eta)*log(1.+eta)-log(zeta)*log(1.+zeta))
		   +ReLi2(-eta)-ReLi2(-zeta)));
    Anu=nu*A;
    vfinite=0.5*A1-nu*A3;
  }
  // t==0
  else {
    // virtual part of the dipole
    Anu = (m02+m12)/(m02-m12)*log(m0/m1);
    vfinite=0.5*(Anu-1.);
  }
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*log(sqr(mgamma)/mprod)+vfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*log(MeV2/mprod)+vfinite);}
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::ReBFF(Energy m1,Energy m2,Energy2 s,double  charge,
			     bool    includegamma,Energy  mgamma) {
  // masses etc
  Energy2 m12(m1*m1),m22(m2*m2),mu(0.5*(s-m12-m22)),mprod(m1*m2);
  // parameters
  double ratio(m1*m2/mu),rho(sqrt((1.-ratio)*(1.+ratio)));
  Energy2 prod(mu*(1.+rho));
  // the finite piece
  double vfinite(mu*rho/s*log(prod/mprod)+0.5*(m12-m22)/s*log(m1/m2)
		 +1./rho*(pi*pi-0.5*log(prod/m12)*log(prod/m22)
			  -0.5*sqr(log((m12+prod)/(m22+prod)))
			  -ReLi2(2.*mu*rho/(m12+prod))
			  -ReLi2(2.*mu*rho/(m22+prod)))-1.);
  // the cut-off piece
  double Anu(log(prod/mprod)/rho),output;
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*log(sqr(mgamma)/mprod)+vfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*log(MeV2/mprod)+vfinite);}
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::BtildeIF(double  beta0   ,double  ombeta0 ,
				double  beta1   ,double  ombeta1 ,
				Energy  en0     ,Energy  en1     ,
				Energy  m0      ,Energy  m1      , 
				Energy2 t       ,double  charge  ,
				Energy  emin    ,bool    includegamma,
				Energy  mgamma) {
  // coefficient of the divergent piece
  Energy2 mprod(m0*m1),nu(0.5*(m0*m0+m1*m1-t));
  double Anu;
  if(nu-mprod>1e-12*GeV2) {
    Energy2 lambda(sqrt((nu-mprod)*(nu+mprod)));
    Anu=nu/lambda*log((lambda+nu)/mprod);
  }
  else
    {Anu=1.;}
  // finite piece
  double rfinite(-0.5*A4single(beta0,ombeta0)-0.5*A4single(beta1,ombeta1)
		 +nu*A4IF(beta0,ombeta0,beta1,ombeta1,en0,en1,m0,m1,t));
  //  assert(isfinite(rfinite));
  // return the answer
  double output; 
  if(includegamma) {
    output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/mgamma)+rfinite);
  }
  else {
    output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/MeV)+rfinite);
  }
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::BtildeFF(double  beta1   ,double  ombeta1 ,
				       double  beta2   ,double  ombeta2 ,
				       Energy  en1     ,Energy  en2     ,
				       Energy  m1      ,Energy  m2      , 
				       Energy2 s       ,double  charge  ,
				       Energy  emin    ,bool    includegamma,
				       Energy  mgamma) {
  // masses etc
  Energy2 m12(m1*m1),m22(m2*m2),mu(0.5*(s-m12-m22)),mprod(m1*m2);
  // parameters
  double ratio(m1*m2/mu),rho(sqrt((1.-ratio)*(1.+ratio)));
  Energy2 prod(mu*(1.+rho));
  // finite piece
  double rfinite(-0.5*A4single(beta1,ombeta1)-0.5*A4single(beta2,ombeta2)
		 +mu*A4FFFull(en1,en2,beta1,beta2,m1,m2,s));
  double Anu(log(prod/mprod)/rho);
  // return the answer
  double output; 
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/mgamma)+rfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/MeV)+rfinite);}
  //  assert(isfinite(output));
  return output;
}

InvEnergy2 YFSFormFactors::A4FFFull(Energy  inen1  ,Energy inen2,
					   double  beta1,double beta2,
					   Energy   inm1  ,Energy inm2,Energy2 s    ) {
  Energy en1(inen1),en2(inen2),m1(inm1),m2(inm2);
  // order the particles so en1>en2
  if(inen1*beta1<inen2*beta2) {
    en1=inen2;
    en2=inen1;
    m1=inm2;
    m2=inm1;
  }
  Energy Delta(en1-en2);
  Energy Omega(en1+en2),delta(m1-m2),omega(m1+m2);
  Energy2 Q2(s-2.*(m1*m1+m2*m2));
  Energy root(sqrt(Delta*Delta+Q2));
  Energy eta[2]={sqrt((en2-m2)*(en2+m2)),sqrt((en1-m1)*(en1+m1))+root};
  if(0.5*(s-m1*m1-m2*m2)>en1*en2){eta[0]=-eta[0];}
  Energy2 root2(sqrt((Q2+omega*omega)*(Q2+delta*delta)));
  double Y[2];
  // various limits
  Energy y[4];
  y[0]=0.5*(root-Omega+(omega*delta+root2)/(root+Delta));
  y[1]=y[0]-root2/(root+Delta);
  y[2]=0.5*(root+Omega+(omega*delta+root2)/(root-Delta));
  y[3]=y[2]-root2/(root-Delta);
  // the Y function at both limits
  for(unsigned int ix=0;ix<2;++ix)
    {Y[ix]=Zij(eta[ix],y[0],y[3])+Zij(eta[ix],y[1],y[0])
	+Zij(eta[ix],y[2],y[1])-Zij(eta[ix],y[2],y[3])
	+0.5*Xijkl(eta[ix],y[0],y[1],y[2],y[3])*Xijkl(eta[ix],y[1],y[2],y[0],y[3]);}
  // the answer
  // the Z function at both limits
  double output(0.);
  if(abs(Delta)>_ecut) {
    output=log(abs((root-Delta)/(root+Delta)))*(+Xijkl(eta[1],y[0],y[3],y[1],y[2])
						-Xijkl(eta[0],y[0],y[3],y[1],y[2]));
  }
  return 1./root2*(output+Y[1]-Y[0]);
}

InvEnergy2 YFSFormFactors::A4IF(double  beta0   ,double  ombeta0 ,
				       double  beta1   ,double  ombeta1 ,
				       Energy  en0  ,Energy en1  , Energy  m0   ,Energy m1   ,
				       Energy2 t) {
  // this is the general function so pick the special case
  if(t>_tcut){
    // rest frame of decaying particle t!=0
    if(abs(en0-m0)<_ecut){return A4IFRest(m0,m1,beta1,ombeta1,en1);}
    // rest frame of decay product t!=0
    else if(abs(en1-m1)<_ecut){return A4IFRest(m1,m0,beta0,ombeta0,en0);}
    // general frame t!=0
    else
      {return A4IFFull(beta0,beta1,en0,en1,m0,m1,t);}
  }
  else {
    // rest frame of decaying particle t=0
    if(abs(en0-m0)<_ecut){return A4IFRestZero(m0,m1);}
    // rest frame of decay products t=0
    else if(abs(en1-m1)<_ecut){return A4IFRestZero(m1,m0);}
    // general frame t=0
    else{return A4IFZero(beta0,beta1,ombeta1,en0,en1,m0,m1);}
  }
}

InvEnergy2 YFSFormFactors::A4IFZero(double  beta0, double beta1, double ombeta1,
				    Energy  en0,
				    Energy en1  , Energy  m0   , Energy m1) {
  Energy  Delta = en0-en1;
  Energy2 mu2 = (m0-m1)*(m0+m1);
  long double z[2]={ beta1*en1/Delta, beta0*en0/Delta-1. };
  long double y[3],xi[3];
  y[0]=en1/Delta;
  y[1]=y[0]-0.5*mu2/sqr(Delta);
  y[2]=-y[0]+2.*m1*m1/mu2;
  for(unsigned int ix = 0; ix < 3; ++ix) {
    if ( ix == 0 ) xi[0]  = -ombeta1*y[0]  / (z[1] - y[0] );
    else           xi[ix] = (z[0] - y[ix]) / (z[1] - y[ix]);
  }
  long double U[2];
  for(unsigned int ix=0;ix<2;++ix) {
    // U[ix] = 0.5*sqr(log(abs((z[ix]-y[0])*(z[ix]-y[1])/(z[ix]-y[2]))))
    //   +log(abs(z[ix]-y[0]))*log(abs(z[ix]-y[0])/sqr(z[ix]-y[1]))
    //   +2.*ReLi2((y[1]-y[0])/(z[ix]-y[0]))
    //   +2.*ReLi2((y[2]-y[1])/(z[ix]-y[1]));
    const long double a = ix==0 ? -ombeta1*y[0] : z[ix]-y[0];
    const long double b = z[ix]-y[1];
    const long double c = z[ix]-y[2];
    const long double A = abs(a*b/c);
    const long double B = abs(a);
    const long double C = B/sqr(b);
    const long double D = (y[1]-y[0])/a;
    const long double E = (y[2]-y[1])/b;
    U[ix] = 0.5*sqr(log(A)) + log(B)*log(C) + 2.*ReLi2(D) + 2.*ReLi2(E);
  }
  return 1./mu2*(log(2.*sqr(Delta)/mu2)*log(abs(xi[1]*xi[2]/xi[0]))+U[1]-U[0]);
}

InvEnergy2 YFSFormFactors::A4IFRest(Energy m0   ,Energy m1, double beta1,
					   double ombeta1, Energy E1) {
  Energy  Mfact0 = m0-E1*ombeta1;
  Energy  Mfact1 = m0-E1*(1.+beta1);
  Energy2 Mfact2 = m0*E1*(1.+beta1)-m1*m1;
  Energy2 Mfact3 = m0*E1*ombeta1-m1*m1;
  Energy2 qprod(m0*E1*beta1);
  return 0.5/qprod*(+log(abs(Mfact0/Mfact1))*log(E1*(1.+beta1)/m0)
		    -2.*log(abs(2.*beta1*E1*Mfact0/m0/m1))*log(E1*(1.+beta1)/m1)
		    +2.*ReLi2(E1/m0*ombeta1)-2.*ReLi2(E1/m0*(1.+beta1))
		    +ReLi2(-0.5*Mfact1/beta1/E1)-ReLi2( 0.5*Mfact0/beta1/E1)
		    +ReLi2( 0.5*Mfact2/qprod   )-ReLi2(-0.5*Mfact3/qprod));
}

InvEnergy2 YFSFormFactors::A4IFFull(Velocity beta0,Velocity beta1,
					   Energy  en0  ,Energy en1  ,
					   Energy  m0   ,Energy m1   , Energy2 t) {
  Energy Delta(en0-en1),Omega(en0+en1),delta(m0-m1),omega(m0+m1);
  Energy  T(sqrt(sqr(Delta)-t)),V(Delta+T);
  Energy2 kappa(sqrt((sqr(omega)-t)*(sqr(delta)-t)));
  long double y[4]={-0.5/T*(T+Omega-(omega*delta+kappa)*V/t),
		    -0.5/T*(T+Omega-(omega*delta-kappa)*V/t),
		    -0.5/T*(T-Omega+(omega*delta+kappa)/V),
		    -0.5/T*(T-Omega+(omega*delta-kappa)/V)};
  long double z[2]={beta1*en1/T,beta0*en0/T-1.};
  double Y[2],lfact(log(abs(V*V/t)));
  for(unsigned int ix=0;ix<2;++ix) {
    Y[ix] = lfact*Xijkl(z[ix],y[0],y[3],y[1],y[2])
      +Zij(z[ix],y[0],y[3])
      +Zij(z[ix],y[1],y[0])
      +Zij(z[ix],y[2],y[1])
      -Zij(z[ix],y[2],y[3])
      +0.5*Xijkl(z[ix],y[0],y[1],y[2],y[3])*Xijkl(z[ix],y[1],y[2],y[0],y[3]);
  }
  return (Y[1]-Y[0])/kappa;
}
