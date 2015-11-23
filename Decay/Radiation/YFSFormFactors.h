// -*- C++ -*-
//
// YFSFormFactors.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_YFSFormFactors_H
#define HERWIG_YFSFormFactors_H
//
// This is the declaration of the YFSFormFactors class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig/Utilities/Maths.h"

namespace Herwig {

using namespace ThePEG;

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
      return -_alpha/Constants::pi*charge*
	((1.+beta1*beta2)/(beta1+beta2)*(+log((1.+beta1)/ombeta1)
					 +log((1.+beta2)/ombeta2)))*log(Emax/Emin);
    else
      return -_alpha/Constants::pi*charge*
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
      return -_alpha/Constants::pi*charge/beta* log((1.+beta)/ombeta)    *log(Emax/Emin);
    else
      return -_alpha/Constants::pi*charge/beta*(log((1.+beta)/ombeta)-2.)*log(Emax/Emin);
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
		      Energy  mgamma=_mgamma);
  
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
		      bool    includegamma=true,Energy  mgamma=_mgamma);
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
			 Energy  mgamma=_mgamma);

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
			 Energy  mgamma=_mgamma);
  //@}
  
  /**
   * Access to the photon mass
   */
  static Energy photonMass() {return _mgamma;}

private:

  /**
   *  Various special cases of the \f$A_4\f$ functions of hep-ph/0302065 for
   *  the initial-final dipole
   */
  //@{
  /**
   * The \f$A_4\f$ function for the full final-final dipole
   * @param inen1 The energy of the first  particle
   * @param inen2 The energy of the second particle
   * @param beta1 Velocity of the first  particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second particle, \f$\beta_2\f$.
   * @param inm1 The mass of the first particle
   * @param inm2 The mass of the second particle
   * @param s   The invariant mass of the charged particles
   */
  static InvEnergy2 A4FFFull(Energy  inen1  ,Energy inen2,
			     double  beta1,double beta2,
			     Energy   inm1  ,Energy inm2,Energy2 s    );
  
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
			 Energy2 t);
  
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
    return -2./mdiff*(sqr(log(m0/m1))+Math::ReLi2(mdiff/sqr(m0)));
  }

  /**
   * The \f$A_4\f$ function for the initial-final dipole with \f$t=0\f$.
   * @param beta0 Velocity of the incoming charged particle, \f$\beta_0\f$
   * @param beta1 Velocity of the outgoing charged particle, \f$\beta_1\f$.
   * @param ombeta1 \f$1-\beta_1\f$ for the outgoing charged particle.
   * @param en0 The energy of the incoming particle
   * @param en1 The energy of the outgoing particle
   * @param m0 The mass of the incoming particle
   * @param m1 The mass of the outgoing particle
   */
  static InvEnergy2 A4IFZero(double  beta0, double beta1, 
			     double ombeta1, Energy  en0,
			     Energy en1  , Energy  m0   , Energy m1);
  
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
			     double ombeta1, Energy E1);
  
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
			     Energy  m0   ,Energy m1   , Energy2 t);
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
  template <typename T> 
  static double Zij(T eta, T yi, T yj) {
    return 2.*Math::ReLi2((yj-yi)/(eta-yi))+0.5*sqr(log(abs((eta-yi)/(eta-yj))));
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
  template <typename T> 
  static double Xijkl(T eta,T yi, T yj, T yk, T yl) {
    return log(abs((eta-yi)*(eta-yj)/(eta-yk)/(eta-yl)));
  }
  //@}

};

}

#endif /* HERWIG_YFSFormFactors_H */
