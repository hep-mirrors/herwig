// -*- C++ -*-
//
// GroupInvariants.h is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#ifndef HERWIG_GroupInvariants_H
#define HERWIG_GroupInvariants_H
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Config/Unitsystem.h"
#include <cassert>

namespace Herwig {
using namespace ThePEG;

namespace GroupInvariants {

  /**
   *  Simple struct for storing the different gauge contributions
   */
  struct GaugeContributions {

    /**
     *  Default Constructor
     */
    GaugeContributions(double inSU3=0.,
		       double inSU2=0., double inU1=0.) 
      : SU3(inSU3),SU2(inSU2),U1(inU1)
    {}
    /**
     * \f$SU(3)\f$
     */
    double SU3;

    /**
     *  \f$SU(2)\f$
     */
    double SU2;

    /**
     * \f$U(1)\f$
     */
    double U1;
  };


  /**
   *  The \f$SU(N)\f$ \f$C_A\f$ 
   */
  inline double C_A(unsigned int N) {
    return N !=1 ? double(N) : 0.;
  }

  /**
   *  The \f$SU(N)\f$ \f$C_F\f$ 
   */
  inline double C_F(unsigned int N) {
    return N !=1 ? 0.5*(double(N*N)-1.)/double(N) : 1.;
  }

  /*
   *  The \f$SU(N)\f$ \f$C_d\f$
   */
  inline double C_d(unsigned int N) {
    return (double(N*N)-4.)/double(N);
  }

  /**
   *  The \f$SU(N)\f$\f$C_1\f$ 
   */
  inline double C_1(unsigned int N) {
    double N2(N*N);
    return 0.25*(N2-1.0)/N2;
  }

  /**
   *  \f$T_F\f$
   */
  inline double T_F(unsigned int N, bool high) {
    if(high) {
      return N !=1 ? 0.5 : 5.0/3.0;
    }
    else {
      return N !=1 ? 0.5 : 20.0/3.0;
    }
  }

  /** 
   *   \f$t_S\f$
   */
  inline double t_S(unsigned int, bool ) {
    return 0.5;
  }

  /**
   * / Number of complex scalars in the fundamental rep. of SU(N)/U(1)
   */
  inline double n_S(unsigned int N, bool high) {
    if(high) {
      if(N==2 || N==1) return 1.0;
      else if(N==3)    return 0.0;
      else assert(false);
    }
    else {
      if(N>=1&&N<=3) return 0.;
      else assert(false);
    }
  }

  /**
   * Number of Dirac Fermions in the fund. rep. of SU(N) (or U(1) for N==1)
   */
  inline double n_F(unsigned int N, bool high) {
    if(high) {
      if(N==1) return 3.0;
      else if(N==2 || N==3) return 6.0;
      else assert(false);
    }
    else {
      if(N==1) return 1.0;
      else if(N==2) return 0.0;
      else if(N==3) return 5.0;
      else assert(false);
    }
  } 
  
  /**
   * Find K_i for gauge group N. high=false for running at mu<mZ
   */
  double K_Factor(unsigned int i,unsigned int N, bool high);

  /**
   *   Find B_i for gauge group N, high energy
   */
  double B_Factor(int i, int N, bool fermion, bool longitudinal);

  /**
   *   Find B_i for gauge group N, low  energy
   */
  double B_Factor_Low(int i, int N, bool fermion, double boostFactor);

  /**
   * Contributions to the Cusps
   */
  GaugeContributions cuspContributions(Energy mu, int K_ORDER, bool high);

  /**
   * Contributions to B, high energy
   */
  GaugeContributions BContributions(Energy mu, int B_ORDER,
				    bool fermion, bool longitudinal);

  /**
   * Contributions to B, low  energy
   */
  GaugeContributions BContributionsLow(Energy mu, int B_ORDER,
				       bool fermion, double boostFactor);

  inline Complex PlusLog(double arg) {
    static const Complex I(0,1.0);
    if (arg>0.0) 
      return log(arg);
    else if (arg<0.0)
      return log(-arg)+I*Constants::pi;
    else 
      assert(false);
  }

  inline Complex MinusLog(double arg) {
    static const Complex I(0,1.0);
    if (arg>0.0) 
      return log(arg);
    else if (arg<0.0)
      return log(-arg)-I*Constants::pi;
    else 
      assert(false);
  }

  /**
   * Number of fermion generations (only used in gauge boson HighCMatching)
   */
  inline double n_g() { return 3.0; }

  /**
   * Number of complex scalars in the fundamental rep. of SU(N)
   */
  inline double nSWeyl(unsigned int N, bool high) {
    if(high) {
      if(N==2 || N==1) return 1.0;
      else if   (N==3) return 0.0;
      else assert(false);
    }
    else {
      if( N==1 || N==3 ) return 0.0;
      else assert(false);
    }
  }

  /**
   * Number of Weyl Fermions in the fundamental rep. of SU(N)
   */
  inline double nFWeyl(unsigned int N, bool high) {
    if(high) {
      if(N==2 || N==3) return 12.0;
      else assert(false);
    }
    else {
      if(N==3)      return 10.0;
      else if(N==1) return  2.0;
      else assert(false);
    }
  }

  inline double TFWeyl(unsigned int) {
    return 0.5;
  }

  inline double tSWeyl(unsigned int) {
    return 0.5;
  }

  inline Complex WFunction(Energy mu, Energy2 s) {
    using Constants::pi;
    assert(abs(s)>ZERO);
    Complex ln = MinusLog(-s/(mu*mu));
    return (-1.0*ln*ln + 3.0*ln+pi*pi/6.0-8.0);
  }

  /**
   * \fX_N\f% function, v is either t or u
   */
  inline Complex XNFunction(unsigned int N, Energy mu, Energy2 s, Energy2 v) {
    using Constants::pi;
    assert(abs(s)>ZERO);
    Complex ls = MinusLog(-s/(mu*mu));
    return (2.0*C_F(N)*WFunction(mu,s) + 
  	    C_A(N)*(2.0*ls*ls - 2.0*MinusLog((s+v)/(mu*mu))*ls - 
  		   11.0/3.0*ls + pi*pi + 85.0/9.0) + 
  	    (2.0/3.0*ls - 10.0/9.0) * TFWeyl(N) * nFWeyl(N,true) + 
  	    (1.0/3.0*ls - 8.0/9.0) * TFWeyl(N) * nSWeyl(N,true));
  }

  /**
   *  \f$\Pi_1\f$ function
   */
  inline Complex PI1_function(Energy mu, Energy2 s) {
    assert(abs(s)>ZERO);
    return ((41.0/6.0)*MinusLog(-s/(mu*mu))-104.0/9.0);
  }

  /**
   *  \f$\tilde{f}\f$ function, v is either t or u
   */
  inline Complex fTildeFunction(Energy mu, Energy2 s, Energy2 v) {
    using Constants::pi;
    assert(abs(s)>ZERO);
    Complex ls  = MinusLog(-s/GeV2), lv = MinusLog(-v/GeV2);
    Complex lsv = MinusLog((s+v)/GeV2);
    return (-2.0*double(s/(s+v))*(lv-ls) + 
   	    double(s*(s+2.0*v)/((s+v)*(s+v))) * ((lv-ls)*(lv-ls) + pi*pi) + 
  	    4.0*MinusLog(-s/(mu*mu))*(lv-lsv));
  }
}
}

#endif // HERWIG_GroupInvariants_H 
