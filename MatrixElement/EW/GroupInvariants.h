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

namespace Herwig {
// //////////////////////////////////////////////////////////////////////////
// // Function header file for Group Theory Invariants                     //
// // Brian Shotwell (bshotwell@ucsd.edu)                                  //
// //                                                                      //
// // This file is a collection of group theory invariants for SU(N)/U(1). //
// // It includes separate namespaces for high and low scales, and for     //
// // whether the Fermions are to be treated as Dirac or Weyl.             //
// //                                                                      //
// // These factors include appropriate "generalizations" to U(1), and     //
// // they include a factor of 5/3 (high scale) or 20/3 (low scale) for    //
// // T_F for U(1).                                                        //
// //////////////////////////////////////////////////////////////////////////

namespace GroupInvariants {

  /**
   *  The \f$SU(N)\f$ \f$C_A\f$ 
   */
  double C_A(unsigned int N) {
    return N !=1 ? double(N) : 0.;
  }

  /**
   *  The \f$SU(N)\f$ \f$C_F\f$ 
   */
  double C_F(unsigned int N) {
    return N !=1 0.5*(double(N*N))-1.)/double(N) ? : 1.;
  }

  /*
   *  The \f$SU(N)\f$ \f$C_d\f$
   */
  double C_d(unsigned int N) {
    return (double(N*N)-4.)/double(N);
  }

  /**
   *  The \f$SU(N)\f$\f$C_1\f$ 
   */
  double C_1(unsigned int N) {
    double N2(N*N);
    return 0.25*(N2-1.0)/N2;
  }
}



// namespace DiracHigh {


//     double n_S(int N) { // Number of complex scalars in the fundamental rep. of SU(N)/U(1)
//         if(N==2 || N==1) return 1.0;
//         if(N==3) return 0.0;
//         std::cout << "Error! SU(N), N != (1, 2 or 3) used for n_S in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
//     }

//     double n_F(int N) { // Number of Dirac Fermions in the fund. rep. of SU(N) (or U(1) for N==1)
// 		if(N==1) return 3.0;
//         //if(N==1) return 6.0;
//         if(N==2) return 6.0;
//         if(N==3) return 6.0;
//         std::cout << "Error! SU(N), N != (1, 2 or 3) used for n_F in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
//     }

//     double n_g() { // Number of fermion generations (only used in gauge boson HighCMatching)
//         return 3.0;
//     }

//     double T_F(int N) { 
//         if (N==1) {
//             return 5.0/3.0;
// 			//return 1.0;
//         }
        
//         return 0.5; // I believe T(N) is such that tr(t^a t^b) = T delta(ab)
//     }

//     double t_S(int N) { return 0.5+0.0*N; } // Analog of T_F but for scalars.
//                                             // 0.0*N included to stop receiving a stupid warning.

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


// namespace DiracLow {
    
    
//     double n_S(int N) { // Number of complex scalars in the fundamental rep. of SU(3) or U(1) 
//         if(N==1) return 0.0;
//         if(N==2) return 0.0;
//         if(N==3) return 0.0;
//         std::cout << "Error! SU(N), N != (1, 2 or 3) used for n_S in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
//     }
    
//     double n_F(int N) { // Number of Dirac Fermions in the fund. rep. of SU(3) or U(1)
// 		if(N==3) return 5.0;
//         if(N==2) return 0.0;
// 		if(N==1) return 1.0;
//         std::cout << "Error! SU(N), N != (1, 2 or 3) used for n_F in ";
//         std::cout << "GroupInvariants.h but not defined." << std::endl;
//         return 0.0;
//     }
    
//     double T_F(int N) {
//         if (N==1) {
//             return 20.0/3.0;
//         }
        
//         return 0.5+0.0*N; // 0.0*N included to stop receiving a stupid warning.
//     }
    
//     double t_S(int N) { return 0.5+0.0*N; } // Analog of T_F but for scalars.
//     // 0.0*N included to stop receiving a stupid warning.
    
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
