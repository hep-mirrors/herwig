// -*- C++ -*-
//
// GroupInvariants.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#include "GroupInvariants.h"
#include "ElectroWeakReweighter.h"

using namespace Herwig;
using namespace GroupInvariants;

double GroupInvariants::K_Factor(unsigned int i,
				 unsigned int N, bool high) {
  using Constants::pi;
  using Constants::zeta3;
  // Find K_i for SU(N) (or, for U(1) if N==1)
  // Relevent for finding the cusp anom. dim.
  if (N!=1 && N!=2 && N!=3) {
    std::cerr << "Error. In AnomDim, function K_Factor, N!={1,2,3}" 
	      << std::endl;
    assert(false);
  }
   
  double CA = C_A(N);
  double CF = C_F(N);
  double nF = n_F(N,high);
  double nS = n_S(N,high);
  double TF = T_F(N,high);
  double tS = t_S(N,high);
  
  if (i==1) {
    return (67.0/9.0-1.0/3.0*pi*pi)*CA - 20.0/9.0*nF*TF - 8.0/9.0*tS*nS;
  }
  else if (i==2) {
    return (245.0/6.0-134.0/27.0*pi*pi+22.0/3.0*zeta3+11.0/45.0*pi*pi*pi*pi)*CA*CA + 
      (-418.0/27.0+40.0/27.0*pi*pi-56.0/3.0*zeta3)*CA*TF*nF + 
      (-55.0/3.0+16.0*zeta3)*CF*TF*nF - 16.0/27.0*TF*TF*nF*nF;
  }
  else {
    std::cerr << "Error. In AnomDim, function K_Factor, i!={1,2}" << std::endl;
    assert(false);
  }
}

GaugeContributions GroupInvariants::cuspContributions(Energy mu, int K_ORDER,
						      bool high) {
  using ThePEG::Constants::pi;
  if (K_ORDER > 3) {
    std::cerr << "Cusp anom dim. requested for K_ORDER>3.\n";
    assert(false);
  }
  // couplings (alphas) for U1, SU2, and SU3
  double a[3];
  if (high) {
    a[0] = ElectroWeakReweighter::coupling()->a1(mu)/(4.0*pi);
    a[1] = ElectroWeakReweighter::coupling()->a2(mu)/(4.0*pi);
    a[2] = ElectroWeakReweighter::coupling()->a3(mu)/(4.0*pi);
  }
  else {
    a[0] = ElectroWeakReweighter::coupling()->aEM(mu)/(4.0*pi);
    a[1] = 0.0;
    a[2] = ElectroWeakReweighter::coupling()->aS(mu)/(4.0*pi);
  }
	
  // gammaN[0] is really 4.0*C_F, but there is a factor of C_F included in the 
  // G1-G3 matrices passed from HighRunning.h to the soft integrator. 
  double gamma1[3];
  gamma1[0] = 4.0;
  gamma1[1] = K_Factor(1,1,high)*gamma1[0];
  gamma1[2] = K_Factor(2,1,high)*gamma1[0];
  double gamma2[3];
  gamma2[0] = 4.0;
  gamma2[1] = K_Factor(1,2,high)*gamma2[0];
  gamma2[2] = K_Factor(2,2,high)*gamma2[0];
  double gamma3[3];
  gamma3[0] = 4.0;
  gamma3[1] = K_Factor(1,3,high)*gamma3[0];
  gamma3[2] = K_Factor(2,3,high)*gamma3[0];
   
  GaugeContributions result;
  if (K_ORDER==0) {
    return result;
  }
  // LO bit
  if (K_ORDER>=1) {
    result.U1 += a[0]*gamma1[0];
    result.SU2 += a[1]*gamma2[0];
    result.SU3 += a[2]*gamma3[0];
  }
  // NLO bit
  if (K_ORDER>=2) {
    result.U1 += a[0]*a[0]*gamma1[1];
    result.SU2 += a[1]*a[1]*gamma2[1];
    result.SU3 += a[2]*a[2]*gamma3[1];
  }
  // NNLO bit
  if (K_ORDER>=3) {
    result.U1 += a[0]*a[0]*a[0]*gamma1[2];
    result.SU2 += a[1]*a[1]*a[1]*gamma2[2];
    result.SU3 += a[2]*a[2]*a[2]*gamma3[2];
  }
  return result;
}

double GroupInvariants::B_Factor(int i, int N, bool fermion,
				 bool longitudinal) {
  using Constants::pi;
  using Constants::zeta3;
  // Find B_i for SU(N) (or, for U(1) if N==1)
  // Relevent for finding the collinear non-cusp anom. dim.
  if (N!=1 && N!=2 && N!=3) {
    std::cerr << "Error. In AnomDim, function B_Factor, N!={1,2,3}\n";
    assert(false);
  }
   
  double CA = C_A(N);
  double CF = C_F(N);
  double nF = n_F(N,true);
  double nS = n_S(N,true);
  double TF = T_F(N,true);
  double tS = t_S(N,true);
  
  if (longitudinal) {
    if (i==1)      return -8.0*CF;
    // Two loop non-cusp not known for scalars
    else if (i==2) return 0.0;
    else assert(false);
  }
  else if (fermion) {
    if (i==1) return -6.0*CF;
    else if (i==2) {
      return (4.0*pi*pi - 48.0*zeta3 - 3.0)*CF*CF + 
	(52.0*zeta3 - 11.0*pi*pi/3.0 - 961.0/27.0)*CA*CF + 
	(4.0*pi*pi/3.0 + 260.0/27.0)*CF*TF*nF + 
	(pi*pi/6.0 + 167.0/54.0)*CF*nS*tS*2.0;
    }
    else
      assert(false);
  }
  else {
    if (i==1) {
      return -2.0*( 11.0*CA/3.0 - 4.0*TF*nF/3.0 - nS*tS/3.0 );
    }
    else if (i==2) {
      return CA*CA*(11.0*pi*pi/9.0+4.0*zeta3-1384.0/27.0) + 
	2.0*CA*nF*TF*(256.0/27.0-2.0*pi*pi/9.0) + 8.0*CF*nF*TF;
    }
    else assert(false);
  }
}

double GroupInvariants::B_Factor_Low(int i, int N, bool fermion,
				     double boostFactor) {
  using Constants::pi;
  using Constants::zeta3;
  // Find B_i for SU(N) (or, for U(1) if N==1)
  // Relevent for finding the collinear non-cusp anom. dim.
  if (N!=1 && N!=2 && N!=3) {
    std::cerr << "Error. In AnomDim, function B_Factor, N!={1,2,3}\n";
    assert(false);
  }
   
  double CA = C_A(N);
  double CF = C_F(N);
  double nF = n_F(N,false);
  double nS = n_S(N,false);
  double TF = T_F(N,false);
  double tS = t_S(N,false);
   
  if (abs(boostFactor)>0.001) {
    if (i==1)      return -4.0*CF;
    // Two loop non-cusp not known for bHQET top, W_L, and W_T fields.
    else if (i==2) return 0.0;
    else assert(false);
  }
  else if (fermion) {
    if (i==1) return -6.0*CF;
    else if (i==2) {
      return (4.0*pi*pi - 48.0*zeta3 - 3.0)*CF*CF + 
	(52.0*zeta3 - 11.0*pi*pi/3.0 - 961.0/27.0)*CA*CF + 
	(4.0*pi*pi/3.0 + 260.0/27.0)*CF*TF*nF + 
	(pi*pi/6.0 + 167.0/54.0)*CF*nS*tS*2.0;
    }
    else
      assert(false);
  }
  // Gluon and Photon are the only things left... use Gauge Boson NonCusps:
  else {
    if (i==1) return -2.0*( 11.0*CA/3.0 - 4.0*TF*nF/3.0 - nS*tS/3.0 );
    else if (i==2) {
      return CA*CA*(11.0*pi*pi/9.0+4.0*zeta3-1384.0/27.0) + 
	2.0*CA*nF*TF*(256.0/27.0-2.0*pi*pi/9.0) + 8.0*CF*nF*TF;
    }
    else 
      assert(false);
  }
}

GaugeContributions GroupInvariants::BContributions(Energy mu,
						   int B_ORDER,
						   bool fermion,
						   bool longitudinal) {
  using Constants::pi;
  // NOTE! THIS RETURNS 2*Gamma+sigma, AND SHOULD BE MULTIPLIED BY 1/2 IN
  // THE COLLINEAR ANOMALOUS DIMENSION INTEGRAND
  if (B_ORDER>2) {
    std::cerr << "Non-cusp collinear anom dim. requested for B_ORDER>2.\n";
    assert(false);
  }
	
  double a[3]; // alpha/(4pi) for U1, SU2, and SU3
	
  a[0] = ElectroWeakReweighter::coupling()->a1(mu)/(4.0*pi);
  a[1] = ElectroWeakReweighter::coupling()->a2(mu)/(4.0*pi);
  a[2] = ElectroWeakReweighter::coupling()->a3(mu)/(4.0*pi);

  GaugeContributions result;
  if (B_ORDER==0) {
    return result;
  }
  if (B_ORDER>=1) {
    result.SU3 += a[2]*B_Factor(1,3,fermion,longitudinal);
    result.SU2 += a[1]*B_Factor(1,2,fermion,longitudinal);
    result.U1  += a[0]*B_Factor(1,1,fermion,longitudinal);
  }
  if (B_ORDER>=2) {
    result.SU3 += a[2]*a[2]*B_Factor(2,3,fermion,longitudinal);
    result.SU2 += a[1]*a[1]*B_Factor(2,2,fermion,longitudinal);
    result.U1  += a[0]*a[0]*B_Factor(2,1,fermion,longitudinal);
  }
  return result;
}

GaugeContributions GroupInvariants::BContributionsLow(Energy mu,
						      int B_ORDER,
						      bool fermion,
						      double boostFactor) {
  using Constants::pi;
  // NOTE! THIS RETURNS 2*Gamma+sigma, AND SHOULD BE MULTIPLIED BY 1/2 IN
  // THE COLLINEAR ANOMALOUS DIMENSION INTEGRAND
  if (B_ORDER>2) {
    std::cerr << "Non-cusp collinear anom dim. requested for B_ORDER>2.\n";
  }
  // alpha/(4pi) for U1, SU2, and SU3
  double a[3]; 
  a[0] = ElectroWeakReweighter::coupling()->aEM(mu)/(4.0*pi);
  a[1] = 0.0;
  a[2] = ElectroWeakReweighter::coupling()->aS(mu)/(4.0*pi);
  GaugeContributions result;
  if (B_ORDER==0) {
    return result;
  }
  if (B_ORDER>=1) {
    result.SU3 += a[2]*B_Factor_Low(1,3,fermion,boostFactor);
    result.SU2 += a[1]*B_Factor_Low(1,2,fermion,boostFactor);
    result.U1  += a[0]*B_Factor_Low(1,1,fermion,boostFactor);
  }
  if (B_ORDER>=2) {
    result.SU3 += a[2]*a[2]*B_Factor_Low(2,3,fermion,boostFactor);
    result.SU2 += a[1]*a[1]*B_Factor_Low(2,2,fermion,boostFactor);
    result.U1  += a[0]*a[0]*B_Factor_Low(2,1,fermion,boostFactor);
  }
  return result;
}
