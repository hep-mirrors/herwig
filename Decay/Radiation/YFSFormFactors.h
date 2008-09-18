// -*- C++ -*-
//
// YFSFormFactors.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_YFSFormFactors_H
#define HERWIG_YFSFormFactors_H
//
// This is the declaration of the YFSFormFactors class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Utilities/Maths.h"

namespace Herwig {

using namespace ThePEG;
using Constants::pi;
using Math::ReLi2;

/**
 * The YFSFormFactors class is a pure static class which implements the various YFS
 * form factors we need for the decays.
 */
class YFSFormFactors {

public:

  /**
   *  The value of \f$\alpha\f$ at \f$q^2=0\f$.
   */
  static const double _alpha;

private:
  /**
   *  The default value of the photon mass
   */
  static const Energy _mgamma;

  /**
   *  The cut-off on the value of \f$t\f$ for the switch to the $t=0$ result
   */
  static const Energy2 _tcut;

  /**
   *  The cut-off on the energy of a particle for it to be considered in its rest frame
   */
  static const Energy _ecut;

public:

  /**
   *  Exponentials of the YFS form factors
   */
  //@{
  /**
   * The exponential of the YFS form factor for the initial-final dipole
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param ombeta0 One minus the velocity of the incoming particle,  \f$1-\beta_0\f$
   * @param ombeta1 One minus the velocity of the outgoing particle,  \f$1-\beta_1\f$
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0  The mass   of the incoming particle
   * @param m1  The mass   of the outgoing particle
   * @param t   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   */
  static double exponentialYFSIF(double  beta0,double  ombeta0,
					double  beta1,double  ombeta1,
					Energy  en0  ,Energy  en1    ,
					Energy  m0   ,Energy  m1     , 
					Energy2 t    ,double  charge ,
					Energy  emin) {
    return exp(YFSIF(beta0,ombeta0,beta1,ombeta1,en0,en1,m0,m1,t,charge,emin));
  }
  /**
   * The exponential of the YFS form factor for the final-final dipole
   * The \f$2\alpha\tilde{B}\f$ function for the final-final dipole
   * @param beta1 Velocity of the first  charged particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second charged particle, \f$\beta_2\f$.
   * @param ombeta1 One minus the velocity of the first  particle,  \f$1-\beta_1\f$
   * @param ombeta2 One minus the velocity of the second particle,  \f$1-\beta_2\f$
   * @param en1 The energy of the first  particle
   * @param en2 The energy of the second particle
   * @param m1  The mass   of the first  particle
   * @param m2  The mass   of the second particle
   * @param s   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   */
  static double exponentialYFSFF(double  beta1, double  ombeta1,
				 double  beta2, double  ombeta2,
				 Energy  en1  , Energy  en2    ,
				 Energy  m1   , Energy  m2     , 
				 Energy2 s    , double  charge ,
				 Energy  emin) {
    return exp(YFSFF(beta1,ombeta1,beta2,ombeta2,en1,en2,m1,m2,s,charge,emin));
  }
  //@}
  /**
   *  The YFS form factors for the initial-final and final-final dipoles
   */
  //@{
  /**
   *  The YFS form factor for the initial-final dipole
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param ombeta0 One minus the velocity of the incoming particle,  \f$1-\beta_0\f$
   * @param ombeta1 One minus the velocity of the outgoing particle,  \f$1-\beta_1\f$
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0  The mass   of the incoming particle
   * @param m1  The mass   of the outgoing particle
   * @param t   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   */
  static double YFSIF(double  beta0   ,double  ombeta0 ,
			     double  beta1   ,double  ombeta1 ,
			     Energy  en0     ,Energy  en1     ,
			     Energy  m0      ,Energy  m1      , 
			     Energy2 t       ,double  charge  ,
			     Energy  emin) {
    return BtildeIF(beta0,ombeta0,beta1,ombeta1,en0,en1,m0,m1,t,charge,emin,false)
      +ReBIF(m0,m1,t,charge,false);
  }
  /**
   *  The YFS form factor for the final-final dipole
   * The \f$2\alpha\tilde{B}\f$ function for the final-final dipole
   * @param beta1 Velocity of the first  charged particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second charged particle, \f$\beta_2\f$.
   * @param ombeta1 One minus the velocity of the first  particle,  \f$1-\beta_1\f$
   * @param ombeta2 One minus the velocity of the second particle,  \f$1-\beta_2\f$
   * @param en1 The energy of the first  particle
   * @param en2 The energy of the second particle
   * @param m1  The mass   of the first  particle
   * @param m2  The mass   of the second particle
   * @param s   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   */
  static double YFSFF(double  beta1   ,double  ombeta1 ,
		      double  beta2   ,double  ombeta2 ,
		      Energy  en1     ,Energy  en2     ,
		      Energy  m1      ,Energy  m2      , 
		      Energy2 s       ,double  charge  ,
		      Energy  emin) {
    return BtildeFF(beta1,ombeta1,beta2,ombeta2,en1,en2,m1,m2,s,charge,emin,false)
      +ReBFF(m1,m2,s,charge,false);
  }
  //@}

  /**
   * Crude average multiplicities for initial-final and final-final dipoles
   */
  //@{
  /** Crude multiplicity for the final-final dipole
   * @param beta1 Velocity of the first  charged particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second charged particle, \f$\beta_2\f$.
   * @param ombeta1 One minus the velocity of the first  particle,  \f$1-\beta_1\f$
   * @param ombeta2 One minus the velocity of the second particle,  \f$1-\beta_2\f$
   * @param charge The product of the charges of the particles in the dipole
   * @param Emin The maximum energy for the integral
   * @param Emax The minimum energy for the integral
   * @param massterms Whether or not to include the mass terms
   */
  static double nbarFF(double beta1, double ombeta1,
			      double beta2, double ombeta2,
			      double charge,
			      Energy Emax , Energy Emin,
			      bool massterms=false) {
    if(!massterms)
      return -_alpha/pi*charge*
	((1.+beta1*beta2)/(beta1+beta2)*(+log((1.+beta1)/ombeta1)
					 +log((1.+beta2)/ombeta2)))*log(Emax/Emin);
    else
      return -_alpha/pi*charge*
	((1.+beta1*beta2)/(beta1+beta2)*(+log((1.+beta1)/ombeta1)
					 +log((1.+beta2)/ombeta2))-2.)*log(Emax/Emin);
  }

  /**
   * Crude multiplicity for the initial-final dipole
   * @param beta Velocity of the outgoing charged particle, \f$\beta\f$
   * @param ombeta One minus the velocity of the outgoing charged  particle,
   *  \f$1-\beta\f$
   * @param charge The product of the charges of the particles in the dipole
   * @param Emin The maximum energy for the integral
   * @param Emax The minimum energy for the integral
   * @param massterms Whether or not to include the mass terms
   */
  //@}
  static double nbarIF(double beta, double ombeta,
		       double charge,
		       Energy Emax , Energy Emin,
		       bool massterms=false) {
    if(!massterms)
      return -_alpha/pi*charge/beta* log((1.+beta)/ombeta)    *log(Emax/Emin);
    else
      return -_alpha/pi*charge/beta*(log((1.+beta)/ombeta)-2.)*log(Emax/Emin);
  }

  /**
   *  The virtual piece for the initial-final and final-final dipoles
   */
  //@{
  /**
   * The \f$2\alpha\mathcal{R}B\f$ function for the initial-final dipole
   * @param m0  The mass   of the incoming particle
   * @param m1  The mass   of the outgoing particle
   * @param t   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param includegamma Include the photon mass terms
   * @param mgamma The photon mass,
   */
  static double ReBIF(Energy  m0      ,Energy  m1      , Energy2 t       ,
		      double  charge  ,bool    includegamma=true,
		      Energy  mgamma=_mgamma) {
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
    return output;
  }
  
  /**
   *  The \f$2\alpha\mathcal{R}B\f$ function for the final-final dipole
   * @param m1  The mass   of the incoming particle
   * @param m2  The mass   of the outgoing particle
   * @param s   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param includegamma Include the photon mass terms
   * @param mgamma The photon mass,
   */
  static double ReBFF(Energy m1,Energy m2,Energy2 s,double  charge,
		      bool    includegamma=true,Energy  mgamma=_mgamma) {
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
    return output;
  }
  //@}

  /**
   *   The real emission terms for initial-final and final-final dipoles
   */
  //@{
  /**
   *  The \f$2\alpha\tilde{B}\f$ function for the initial-final dipole
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param ombeta0 One minus the velocity of the incoming particle,  \f$1-\beta_0\f$
   * @param ombeta1 One minus the velocity of the outgoing particle,  \f$1-\beta_1\f$
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0  The mass   of the incoming particle
   * @param m1  The mass   of the outgoing particle
   * @param t   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   * @param includegamma Include the photon mass terms
   * @param mgamma The photon mass,
   */
  static double BtildeIF(double  beta0   ,double  ombeta0 ,
			 double  beta1   ,double  ombeta1 ,
			 Energy  en0     ,Energy  en1     ,
			 Energy  m0      ,Energy  m1      , 
			 Energy2 t       ,double  charge  ,
			 Energy  emin    ,bool    includegamma=true,
			 Energy  mgamma=_mgamma) {
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
    // return the answer
    double output; 
    if(includegamma) {
      output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/mgamma)+rfinite);
    }
    else {
      output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/MeV)+rfinite);
    }
    return output;
  }
  
  /**
   * The \f$2\alpha\tilde{B}\f$ function for the final-final dipole
   * @param beta1 Velocity of the first  charged particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second charged particle, \f$\beta_2\f$.
   * @param ombeta1 One minus the velocity of the first  particle,  \f$1-\beta_1\f$
   * @param ombeta2 One minus the velocity of the second particle,  \f$1-\beta_2\f$
   * @param en1 The energy of the first  particle
   * @param en2 The energy of the second particle
   * @param m1  The mass   of the first  particle
   * @param m2  The mass   of the second particle
   * @param s   The invariant mass of the charged particles
   * @param charge The product of the charges of the particles in the dipole
   * @param emin The minimum photon energy
   * @param includegamma Include the photon mass terms
   * @param mgamma The photon mass,
   */
  static double BtildeFF(double  beta1   ,double  ombeta1 ,
			 double  beta2   ,double  ombeta2 ,
			 Energy  en1     ,Energy  en2     ,
			 Energy  m1      ,Energy  m2      , 
			 Energy2 s       ,double  charge  ,
			 Energy  emin    ,bool    includegamma=true,
			 Energy  mgamma=_mgamma) {
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
    return output;
  }
  //@}
  
private:

  /**
   *  Various special cases of the \f$A_4\f$ functions of hep-ph/0302065 for
   *  the initial-final dipole
   */
  //@{
  /**
   * The \f$A_4\f$ function for the full final-final dipole
   * @param en1 The energy of the first  particle
   * @param en2 The energy of the second particle
   * @param beta1 Velocity of the first  particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second particle, \f$\beta_2\f$.
   * @param m1 The mass of the first particle
   * @param m2 The mass of the second particle
   * @param s   The invariant mass of the charged particles
   */
  static InvEnergy2 A4FFFull(Energy  inen1  ,Energy inen2,
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
  
  /**
   *  The \f$A_4\f$ function of hep-ph0302065 using the special cases where
   *  necessary for numerical stability
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param ombeta0 One minus the velocity of the incoming particle,  \f$1-\beta_0\f$
   * @param ombeta1 One minus the velocity of the outgoing particle,  \f$1-\beta_1\f$
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   * @param t  The invariant mass of the charged particles
   */
  static InvEnergy2 A4IF(double  beta0   ,double  ombeta0 ,
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
      else{return A4IFZero(beta0,beta1,en0,en1,m0,m1);}
    }
  }
  
  /**
   * The \f$A_4\f$ function of hep-ph/0302065 for a single particle, without the \f$1/p^2\f$
   * pre-factor
   */
  static double A4single(double beta,double ombeta) {
    if(beta>0.01) return log(ombeta/(1.+beta))/beta;
    else          return -2.-2./3.*sqr(beta)*(1+0.6*sqr(beta));
  }

  /**
   * The \f$A_4\f$ function of hep-ph/0302065 for the initial-final dipole with \f$t=0\f$ in the
   * rest frame.
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   */
  static InvEnergy2 A4IFRestZero(Energy  m0, Energy m1) {
    Energy2 mdiff(m0*m0-m1*m1);
    return -2./mdiff*(sqr(log(m0/m1))+ReLi2(mdiff/sqr(m0)));
  }

  /**
   * The \f$A_4\f$ function for the initial-final dipole with \f$t=0\f$.
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   */
  static InvEnergy2 A4IFZero(double  beta0, double beta1, Energy  en0,
			     Energy en1  , Energy  m0   , Energy m1) {
    Energy  Delta(en0-en1);
    Energy2 mu2((m0-m1)*(m0+m1));
    long double z[2]={ beta1*en1/Delta, beta0*en0/Delta-1. };
    long double y[3],xi[3];
    y[0]=en1/Delta;
    y[1]=y[0]-0.5*mu2/sqr(Delta);
    y[2]=-y[0]+2.*m1*m1/mu2;
    for(unsigned int ix = 0; ix < 3; ++ix) {
      xi[ix] = (z[0] - y[ix]) / (z[1] - y[ix]);
    }
    long double U[2];
    for(unsigned int ix=0;ix<2;++ix)
      {U[ix] = 0.5*sqr(log(abs((z[ix]-y[0])*(z[ix]-y[1])/(z[ix]-y[2]))))
	  +log(abs(z[ix]-y[0]))*log(abs(z[ix]-y[0])/sqr(z[ix]-y[1]))
	  +2.*ReLi2((y[1]-y[0])/(z[ix]-y[0]))
	  +2.*ReLi2((y[2]-y[1])/(z[ix]-y[1]));}
    return 1./mu2*(log(2.*sqr(Delta)/mu2)*log(abs(xi[1]*xi[2]/xi[0]))+U[1]-U[0]);
  }
  
  /**
   * The \f$A_4\f$ function for the initial-final dipole in the rest frame of 
   * the decaying particle
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   * @param beta1 The velocity of the decay product
   * @param ombeta1 \f$1-\beta\f$ for the decay product
   * @param E1 The energy of the outgoing particle
   */
  static InvEnergy2 A4IFRest(Energy m0   ,Energy m1, double beta1,
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
  
  /**
   * The \f$A_4\f$ function for the full initial-final dipole
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   * @param t  The invariant mass of the charged particles
   */
  static InvEnergy2 A4IFFull(Velocity beta0,Velocity beta1,
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
  //@}

  /**
   *  Functions from hep-ph/9606429 for the calculation of the \f$A_4\f$ functions
   */
  //@{
  /**
   * The function \f$Z_{ij}(\eta)\f$ from hep-ph/9606429 for the evaluation of the \f$A_4\f$ 
   * function.
   * @param eta The value of \f$\eta\f$
   * @param yi The value of \f$y_i\f$ 
   * @param yj The value of \f$y_j\f$ 
   */
  template <typename T> static double Zij(T eta, T yi, T yj) {
    return 2.*ReLi2((yj-yi)/(eta-yi))+0.5*sqr(log(abs((eta-yi)/(eta-yj))));
  }
  
  /**
   * The function \f$X^{ij}_{kl}(\eta)\f$  from hep-ph/9606429 for the evaluation of the \f$A_4\f$ 
   * function.
   * @param eta The value of \f$\eta\f$
   * @param yi The value of \f$y_i\f$ 
   * @param yj The value of \f$y_j\f$ 
   * @param yk The value of \f$y_k\f$ 
   * @param yl The value of \f$y_l\f$ 
   */
  template <typename T> static double Xijkl(T eta,T yi, T yj,
					    T yk, T yl) {
    return log(abs((eta-yi)*(eta-yj)/(eta-yk)/(eta-yl)));
  }
  //@}

private:

  /** @name Standard constructors and assignment are private as this is a pure static 
   * class
   */
  //@{
  /**
   * The default constructor.
   */
  YFSFormFactors();

  /**
   * The copy constructor.
   */
  YFSFormFactors(const YFSFormFactors &);

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  YFSFormFactors & operator=(const YFSFormFactors &);
  //@}

};

}

#endif /* HERWIG_YFSFormFactors_H */
