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
}



// namespace DiracHigh {


//     double n_g() { // Number of fermion generations (only used in gauge boson HighCMatching)
//         return 3.0;
//     }


// }


// namespace WeylHigh {
	
// 	double n_S(int N) { // Number of complex scalars in the fundamental rep. of SU(N)/U(1)
// 		if(N==2 || N==1) return 1.0;
// 		if(N==3) return 0.0;
// 		std::cout << "Error! SU(N), N != (1, 2 or 3) used for n_S in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
// 		return 0.0;
// 	}
// 	double n_F(int N) { // Number of Weyl Fermions in the fundamental rep. of SU(N)
// 		if(N==2) return 12.0;
// 		if(N==3) return 12.0;
// 		std::cout << "Error! SU(N), N != (2 or 3) used for n_F in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
// 		return 0.0;
// 	}
//     double n_g() { // Number of fermion generations (only used in gauge boson HighCMatching)
//         return 3.0;
//     }
// 	double T_F(int N) { return 0.5+0.0*N; } // I believe T(N) is such that tr(t^a t^b) = T delta(ab)
// 	// 0.0*N included to stop receiving a stupid warning.
	
// 	double t_S(int N) { return 0.5+0.0*N; } // Analog of T_F but for scalars.
// 	// 0.0*N included to stop receiving a stupid warning.
	
// }

// namespace WeylLow {
	
// 	double n_S(int N) { // Number of complex scalars in the fundamental rep. of SU(N)
//         if(N==1) return 0.0;
//         if(N==3) return 0.0;
//         std::cout << "Error! SU(N), N != (1 or 3) used for n_S in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
// 	}
// 	double n_F(int N) { // Number of Weyl Fermions in the fundamental rep. of SU(N)
// 		if(N==3) return 10.0;
// 		if(N==1) return 2.0;
//         std::cout << "Error! SU(N), N != (1 or 3) used for n_F in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
// 	}
// 	double T_F(int N) { return 0.5+0.0*N; } // I believe T(N) is such that tr(t^a t^b) = T delta(ab)
// 	// 0.0*N included to stop receiving a stupid warning.
	
// 	double t_S(int N) { return 0.5+0.0*N; } // Analog of T_F but for scalars.
// 	// 0.0*N included to stop receiving a stupid warning.
	
// }


// #endif // GROUP_INVARIANTS_H



}

#endif // HERWIG_GroupInvariants_H 
