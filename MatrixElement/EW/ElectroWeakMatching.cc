// -*- C++ -*-
//
// ElectroWeakMatching.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#include "ElectroWeakMatching.h"
#include "ElectroWeakReweighter.h"
#include "GroupInvariants.h"
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace Herwig;
using namespace ElectroWeakMatching;
using namespace GroupInvariants;
using namespace EWProcess;

boost::numeric::ublas::matrix<Complex>
ElectroWeakMatching::electroWeakMatching(Energy mu,
					 Energy2 s, Energy2 t, Energy2 u,
					 Herwig::EWProcess::Process process,
					 bool oneLoop) {
  static const Complex I(0,1.0);
  using Constants::pi;
  Complex T = getT(s,t);
  Complex U = getU(s,u);
   // Z-Couplings
   double g_Lu = ElectroWeakReweighter::coupling()->g_Lu(mu);
   double g_Ld = ElectroWeakReweighter::coupling()->g_Ld(mu);
   double g_Le = ElectroWeakReweighter::coupling()->g_Le(mu);
   double g_Lnu = ElectroWeakReweighter::coupling()->g_Lnu(mu);
   double g_Ru = ElectroWeakReweighter::coupling()->g_Ru(mu);
   double g_Rd = ElectroWeakReweighter::coupling()->g_Rd(mu);
   double g_Re = ElectroWeakReweighter::coupling()->g_Re(mu);
   double g_W = ElectroWeakReweighter::coupling()->g_W(mu);
   double g_phiPlus = ElectroWeakReweighter::coupling()->g_phiPlus(mu);
   // Weinberg Angle:
   double cos2 = ElectroWeakReweighter::coupling()->Cos2thW(mu);
   double sin2 = 1.0-cos2;
   double cos = sqrt(cos2);
   double sin = sqrt(sin2);

   boost::numeric::ublas::matrix<Complex> R0,G2,Dw,Dz;
   
   switch (process) {
   case QQQQ:
   case QQQQiden:
   case QtQtQQ:
     {
       unsigned int numGauge = 4, numBrokenGauge = 12;
       boost::numeric::ublas::matrix<Complex> R0=boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       boost::numeric::ublas::matrix<Complex> G2=boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       boost::numeric::ublas::matrix<Complex> Dw=boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       boost::numeric::ublas::matrix<Complex> Dz=boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,1) = R0(1,1) = R0(2,1) = R0(3,1) = 1.0;
       R0(0,0) = R0(3,0) = 0.25;
       R0(1,0) = R0(2,0) = -0.25;
       R0(4,0) = R0(5,0) = 0.5;
       R0(6,3) = R0(7,3) = R0(8,3) = R0(9,3) = 1.0;
       R0(6,2) = R0(9,2) = 0.25;
       R0(7,2) = R0(8,2) = -0.25;
       R0(10,2) = R0(11,2) = 0.5;
       if (oneLoop) {
	 double g11 = g_Lu;
	 double g12 = g_Ld;
	 double g21 = g_Lu;
	 double g22 = g_Ld;
	 for(unsigned int ix=0;ix<numBrokenGauge;++ix) {
	   Dw(ix,ix) = 0.5*I*pi;
	 }
	 Complex w1 = -0.5*(T-U);
	 Complex w2 = -0.5*(T+U);
	 for(unsigned int ix=0;ix<numBrokenGauge;ix+=6) {
	   Dw(ix+0,ix+0) += w1;
	   Dw(ix+3,ix+3) += w1;
	   Dw(ix+1,ix+1) +=-w1;
	   Dw(ix+2,ix+2) +=-w1;
	   Dw(ix+4,ix+4) += w2; 
	   Dw(ix+5,ix+5) += w2;
	 }
	 Complex z1 = 2.0*g11*g21*(T-U) - I*pi*(g11*g11+g21*g21);
	 Complex z2 = 2.0*g21*g12*(T-U) - I*pi*(g21*g21+g12*g12);
	 Complex z3 = 2.0*g22*g11*(T-U) - I*pi*(g22*g22+g11*g11);
	 Complex z4 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Complex z5 = (g11*g21+g12*g22)*T - (g21*g12+g11*g22)*U 
	   - 0.5*I*pi*(g11*g11+g12*g12+g21*g21+g22*g22);
	 for(unsigned int ix=0;ix<numBrokenGauge;ix+=6) {
	   Dz(ix+0,ix+0) = z1;
	   Dz(ix+1,ix+1) = z2;
	   Dz(ix+2,ix+2) = z3;
	   Dz(ix+3,ix+3) = z4;
	   Dz(ix+4,ix+4) = z5;
	   Dz(ix+5,ix+5) = z5;
	 }  
	 boost::numeric::ublas::matrix<Complex> gam2 = Gamma2(U,T);
	 G2(0,0) += gam2(0,0);
	 G2(0,1) += gam2(0,1);
	 G2(1,0) += gam2(1,0);
	 G2(1,1) += gam2(1,1);
	 G2(2,2) += gam2(0,0);
	 G2(2,3) += gam2(0,1);
	 G2(3,2) += gam2(1,0);
	 G2(3,3) += gam2(1,1);
       }
     }
     break;
   case QQUU:
   case QtQtUU:
   case QQtRtR:
     {
       unsigned int numGauge = 2, numBrokenGauge = 4;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = R0(2,1) = R0(3,1) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lu;
	 double g12 = g_Ld;
	 //double g21 = g_Ru;
	 double g22 = g_Ru;
         
	 Complex w1 = 0.25*I*pi;
	 for(unsigned int ix=0;ix<numBrokenGauge;++ix) Dw(ix,ix) = w1;
	 
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 for(unsigned int ix=0;ix<numBrokenGauge;ix+=2) {
	   Dz(ix+0,ix+0) = z1;
	   Dz(ix+1,ix+1) = z2;
	 }
	 G2 = Gamma2Singlet();
       }
     }
     break;
   case QQDD:
   case QtQtDD:
     {
       unsigned int numGauge = 2, numBrokenGauge = 4;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = R0(2,1) = R0(3,1) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lu;
	 double g12 = g_Ld;
	 //double g21 = g_Rd;
	 double g22 = g_Rd;
         
	 Complex w1 = 0.25*I*pi;
	 for(unsigned int ix=0;ix<numBrokenGauge;++ix) Dw(ix,ix) = w1;
         
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 for(unsigned int ix=0;ix<numBrokenGauge;ix+=2) {
	   Dz(ix+0,ix+0) = z1;
	   Dz(ix+1,ix+1) = z2;
	 }
	 G2 = Gamma2Singlet();
       }
     }
     break;
   case QQLL:
     {
       unsigned int numGauge = 2, numBrokenGauge = 6;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,1) = R0(1,1) = R0(2,1) = R0(3,1) = 1.0;
       R0(0,0) = R0(3,0) = 0.25;
       R0(1,0) = R0(2,0) = -0.25;
       R0(4,0) = R0(5,0) = 0.5;
       if (oneLoop) {
	 double g11 = g_Lu;
	 double g12 = g_Ld;
	 double g21 = g_Lnu;
	 double g22 = g_Le;
	 
	 for (unsigned int i=0; i<6; ++i) {
	   Dw(i,i) = 0.5*I*pi;
	 }
	 Complex w1 = (-1.0/2.0)*(T-U);
	 Complex w2 = (-1.0/2.0)*(T+U);
	 Dw(0,0) += w1;
	 Dw(3,3) += w1;
	 Dw(1,1) += -1.0*w1;
	 Dw(2,2) += -1.0*w1;
	 Dw(4,4) += w2; 
	 Dw(5,5) += w2;
         
	 Complex z1 = 2.0*g11*g21*(T-U) - I*pi*(g11*g11+g21*g21);
	 Complex z2 = 2.0*g21*g12*(T-U) - I*pi*(g21*g21+g12*g12);
	 Complex z3 = 2.0*g22*g11*(T-U) - I*pi*(g22*g22+g11*g11);
	 Complex z4 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Complex z5 = (g11*g21+g12*g22)*T - (g21*g12+g11*g22)*U 
	   - 0.5*I*pi*(g11*g11+g12*g12+g21*g21+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
	 Dz(2,2) = z3;
	 Dz(3,3) = z4;
	 Dz(4,4) = Dz(5,5) = z5;
         
	 G2 = Gamma2(U,T);
       }
     }
     break;
   case QQEE:
     {
       unsigned int numGauge = 1, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lu;
	 double g12 = g_Ld;
	 //double g21 = g_Re;
	 double g22 = g_Re;
            
	 Complex w1 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = w1;
            
         Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
            
	 G2(0,0) = Gamma2Singlet()(0,0);
       }
     }
     break;
   case UUUU:
   case UUUUiden:
   case tRtRUU:
     {
       unsigned int numGauge = 2, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,1) = 1.0;
       if (oneLoop) {
	 double g11 = g_Ru;
	 //double g12 = g_Ru;
	 //double g21 = g_Ru;
	 double g22 = g_Ru;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = Dz(1,1) = z1;
       }
     }
     break;
   case UUDD:
   case tRtRDD:
     {
       unsigned int numGauge = 2, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,1) = 1.0;
       if (oneLoop) {
	 double g11 = g_Ru;
	 //double g12 = g_Ru;
	 //double g21 = g_Rd;
	 double g22 = g_Rd;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = Dz(1,1) = z1;
       }
     }
     break;
   case UULL:
     {
       unsigned int numGauge = 1, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lnu;
	 double g12 = g_Le;
	 //double g21 = g_Ru;
	 double g22 = g_Ru;
         
	 Complex w1 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = w1;
            
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
            
	 G2(0,0) = Gamma2Singlet()(0,0);
       }
     }
     break;
   case UUEE:
     {
       unsigned int numGauge = 1, numBrokenGauge = 1;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Ru;
	 //double g12 = g_Ru;
	 //double g21 = g_Re;
	 double g22 = g_Re;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = z1;
       }
     }
     break;
   case DDDD:
   case DDDDiden:
     {
       unsigned int numGauge = 2, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,1) = 1.0;
       if (oneLoop) {
	 double g11 = g_Rd;
	 //double g12 = g_Rd;
	 //double g21 = g_Rd;
	 double g22 = g_Rd;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = Dz(1,1) = z1;
       }
     }
     break;
   case DDLL:
     {
       unsigned int numGauge = 1, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lnu;
	 double g12 = g_Le;
	 //double g21 = g_Rd;
	 double g22 = g_Rd;
         
	 Complex w1 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = w1;
         
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
         
	 G2(0,0) = Gamma2Singlet()(0,0);
       }
     }
     break;
   case DDEE:
     {
       unsigned int numGauge = 1, numBrokenGauge = 1;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0 *= 0.0; Dw = Dz *= 0.0;
       R0(0,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Rd;
	 //double g12 = g_Rd;
	 //double g21 = g_Re;
	 double g22 = g_Re;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = z1;
       }
     }
     break;
   case LLLL:
   case LLLLiden:
     {
       unsigned int numGauge = 2, numBrokenGauge = 6;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,1) = R0(1,1) = R0(2,1) = R0(3,1) = 1.0;
       R0(0,0) = R0(3,0) = 0.25;
       R0(1,0) = R0(2,0) = -0.25;
       R0(4,0) = R0(5,0) = 0.5;
       if (oneLoop) {
	 double g11 = g_Lnu;
	 double g12 = g_Le;
	 double g21 = g_Lnu;
	 double g22 = g_Le;
         
	 for (int i=0; i<6; i++) {
	   Dw(i,i) = 0.5*I*pi;
	 }
	 Complex w1 = (-1.0/2.0)*(T-U);
	 Complex w2 = (-1.0/2.0)*(T+U);
	 Dw(0,0) += w1;
	 Dw(3,3) += w1;
	 Dw(1,1) += -1.0*w1;
	 Dw(2,2) += -1.0*w1;
	 Dw(4,4) += w2; 
	 Dw(5,5) += w2;
         
	 Complex z1 = 2.0*g11*g21*(T-U) - I*pi*(g11*g11+g21*g21);
	 Complex z2 = 2.0*g21*g12*(T-U) - I*pi*(g21*g21+g12*g12);
	 Complex z3 = 2.0*g22*g11*(T-U) - I*pi*(g22*g22+g11*g11);
	 Complex z4 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Complex z5 = (g11*g21+g12*g22)*T - (g21*g12+g11*g22)*U 
	   - 0.5*I*pi*(g11*g11+g12*g12+g21*g21+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
	 Dz(2,2) = z3;
	 Dz(3,3) = z4;
	 Dz(4,4) = Dz(5,5) = z5;
         
	 G2 = Gamma2(U,T);
       }
     }
     break;
   case LLEE:
     {
       unsigned int numGauge = 1, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(1,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Lnu;
	 double g12 = g_Le;
	 //double g21 = g_Re;
	 double g22 = g_Re;
	 
	 Complex w1 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = w1;
         
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Complex z2 = 2.0*g12*g22*(T-U) - I*pi*(g12*g12+g22*g22);
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
         
	 G2(0,0) = Gamma2Singlet()(0,0);
       }
     }
     break;
   case EEEE:
   case EEEEiden:
     {
       unsigned int numGauge = 1, numBrokenGauge = 1;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 1.0;
       if (oneLoop) {
	 double g11 = g_Re;
	 //double g12 = g_Re;
	 //double g21 = g_Re;
	 double g22 = g_Re;
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = 2.0*g11*g22*(T-U) - I*pi*(g11*g11+g22*g22);
	 Dz(0,0) = z1;
       }
     }
     break;
   case QQWW:
   case LLWW:
     {
       unsigned int numGauge = 5, numBrokenGauge = 20;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 1.0; R0(0,1) = 0.5;
       R0(1,0) = 1.0; R0(1,1) = -0.5;
       R0(2,0) = cos2; R0(2,2) = -0.5*sin*cos; R0(2,3) = -0.5*sin*cos; R0(2,4) = sin2;
       R0(3,0) = sin*cos; R0(3,2) = 0.5*cos2; R0(3,3) = -0.5*sin2; R0(3,4) = -sin*cos;
       R0(4,0) = sin*cos; R0(4,2) = -0.5*sin2; R0(4,3) = 0.5*cos2; R0(4,4) = -sin*cos;
       R0(5,0) = sin2; R0(5,2) = 0.5*sin*cos; R0(5,3) = 0.5*sin*cos; R0(5,4) = cos2;
       R0(6,0) = 1.0; R0(6,1) = -0.5;
       R0(7,0) = 1.0; R0(7,1) = 0.5;
       R0(8,0) = cos2; R0(8,2) = 0.5*sin*cos; R0(8,3) = 0.5*sin*cos; R0(8,4) = sin2;
       R0(9,0) = sin*cos; R0(9,2) = -0.5*cos2; R0(9,3) = 0.5*sin2; R0(9,4) = -sin*cos;
       R0(10,0) = sin*cos; R0(10,2) = 0.5*sin2; R0(10,3) = -0.5*cos2; R0(10,4) = -sin*cos;
       R0(11,0) = sin2; R0(11,2) = -0.5*sin*cos; R0(11,3) = -0.5*sin*cos; R0(11,4) = cos2;
       R0(12,1) = -cos/sqrt(2.0); R0(12,3) = -sin/sqrt(2.0);
       R0(13,1) = -sin/sqrt(2.0); R0(13,3) = cos/sqrt(2.0);
       R0(14,1) = cos/sqrt(2.0); R0(14,2) = -sin/sqrt(2.0);
       R0(15,1) = sin/sqrt(2.0); R0(15,2) = cos/sqrt(2.0);
       R0(16,1) = -cos/sqrt(2.0); R0(16,2) = -sin/sqrt(2.0);
       R0(17,1) = -sin/sqrt(2.0); R0(17,2) = cos/sqrt(2.0);
       R0(18,1) = cos/sqrt(2.0); R0(18,3) = -sin/sqrt(2.0);
       R0(19,1) = sin/sqrt(2.0); R0(19,3) = cos/sqrt(2.0);
       if (oneLoop) {
	 double gW = g_W;
	 double g1(0.),g2(0.);
	 if (process==QQWW) {
	   g1 = g_Lu;
	   g2 = g_Ld;
	 }
	 else if (process==LLWW) {
	   g1 = g_Lnu;
	   g2 = g_Le;
	 }
            
	 Complex w1 = T-U+5.0/4.0*I*pi;
	 Complex w2 = -T+U+5.0/4.0*I*pi;
	 Complex w3 = -0.5*(T+U) + 3.0/4.0*I*pi;
	 Complex w4 = 0.25*I*pi;
	 Dw(0,0) = Dw(7,7) = w1;
	 Dw(1,1) = Dw(6,6) = w2;
	 for (unsigned int i=12; i<20; i++) {
	   Dw(i,i) = w3;
	 }
	 Dw(2,2) = Dw(3,3) = Dw(4,4) = Dw(5,5) = w4;
	 Dw(8,8) = Dw(9,9) = Dw(10,10) = Dw(11,11) = w4;
         
	 Complex z1 = 2.0*g1*gW*(U-T) - I*pi*(g1*g1+gW*gW);
	 Complex z2 = 2.0*g1*gW*(T-U) - I*pi*(g1*g1+gW*gW);
	 Complex z3 = 2.0*g2*gW*(U-T) - I*pi*(g2*g2+gW*gW);
	 Complex z4 = 2.0*g2*gW*(T-U) - I*pi*(g2*g2+gW*gW);
	 Complex z5 = -(g2*gW)*T + (g1*gW)*U - I*pi*(g1*g2+g1*gW-g2*gW);
	 Complex z6 = (g1*gW)*T - (g2*gW)*U - I*pi*(g1*g2+g1*gW-g2*gW);
	 Complex z7 = -I*pi*g1*g1;
	 Complex z8 = -I*pi*g2*g2;
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
	 Dz(2,2) = Dz(3,3) = Dz(4,4) = Dz(5,5) = z7;
	 Dz(6,6) = z3;
	 Dz(7,7) = z4;
	 Dz(8,8) = Dz(9,9) = Dz(10,10) = Dz(11,11) = z8;
	 Dz(12,12) = Dz(13,13) = Dz(16,16) = Dz(17,17) = z5;
	 Dz(14,14) = Dz(15,15) = Dz(18,18) = Dz(19,19) = z6;
         
	 G2 = Gamma2w(U,T);
       }
     }
     break;
   case QQPhiPhi:
   case LLPhiPhi:
     {
       unsigned int numGauge = 2, numBrokenGauge = 14;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 0.25;		R0(0,1) = 1.0;
       R0(1,0) = -1.0/8.0;	R0(1,1) = 0.5;
       R0(2,0) = I/8.0;		R0(2,1) = -I/2.0;
       R0(3,0) = -I/8.0;    R0(3,1) = I/2.0;
       R0(4,0) = -1.0/8.0;	R0(4,1) = 1.0/2.0;
       R0(5,0) = -1.0/4.0;	R0(5,1) = 1.0;
       R0(6,0) = 1.0/8.0;	R0(6,1) = 1.0/2.0;
       R0(7,0) = -I/8.0;    R0(7,1) = -I/2.0;
       R0(8,0) = I/8.0;		R0(8,1) = I/2.0;
       R0(9,0) = 1.0/8.0; 	R0(9,1) = 1.0/2.0;
       R0(10,0) = -1.0/(2.0*sqrt(2.0));
       R0(11,0) = I/(2.0*sqrt(2.0));
       R0(12,0) = -1.0/(2.0*sqrt(2.0));
       R0(13,0) = -I/(2.0*sqrt(2.0));
       
       if (oneLoop) {
	 double g1(0.),g2(0.);
	 if (process==QQWW) {
	   g1 = g_Lu;
	   g2 = g_Ld;
	 }
	 else if (process==LLWW) {
	   g1 = g_Lnu;
	   g2 = g_Le;
	 }
	 else
	   assert(false);
	 double g3 = g_phiPlus;
         
	 Complex w0 = 0.25*I*pi;
	 Complex w1 = 0.5*(T-U) + 0.5*I*pi;
	 Complex w2 = -0.5*(T-U) + 0.5*I*pi;
	 Complex w3 = 0.25*I*(T-U);
	 Complex w4 = -0.25*(T+U) + 0.25*I*pi;
	 Dw(0,0) = w2;
	 Dw(1,1) = w0; Dw(1,2) = w3; Dw(1,3) = -w3; Dw(1,4) = w0;
	 Dw(2,1) = -w3; Dw(2,2) = w0; Dw(2,3) = -w0; Dw(2,4) = -w3;
	 Dw(3,1) = w3; Dw(3,2) = -w0; Dw(3,3) = w0; Dw(3,4) = w3;
	 Dw(4,1) = w0; Dw(4,2) = w3; Dw(4,3) = -w3; Dw(4,4) = w0;
	 Dw(5,5) = w1;
	 Dw(6,6) = w0; Dw(6,7) = -w3; Dw(6,8) = w3; Dw(6,9) = w0;
	 Dw(7,6) = w3; Dw(7,7) = w0; Dw(7,8) = -w0; Dw(7,9) = w3;
	 Dw(8,6) = -w3; Dw(8,7) = -w0; Dw(8,8) = w0; Dw(8,9) = -w3;
	 Dw(9,6) = w0; Dw(9,7) = -w3; Dw(9,8) = w3; Dw(9,9) = w0;
	 Dw(10,10) = w4; Dw(10,11) = I*w4;
	 Dw(11,10) = -I*w4; Dw(11,11) = w4;
	 Dw(12,12) = w4; Dw(12,13) = -I*w4;
	 Dw(13,12) = I*w4; Dw(13,13) = w4;
         
	 Complex z1 = 2.0*g3*g1*(T-U) - I*pi*(g3*g3+g1*g1);
	 Complex z2 = 2.0*g3*g2*(T-U) - I*pi*(g3*g3+g2*g2);
	 Complex z3 = -I*pi*g1*g1;
	 Complex z4 = 0.5*I*g1*(T-U);
	 Complex z5 = 0.25*I*pi;
	 Complex z6 = -I*pi*g2*g2;
	 Complex z7 = 0.5*I*g2*(T-U);
	 Complex z8 = g3*g1*T-g3*g2*U-I*pi*(g1*g2-g2*g3+g1*g3);
	 Complex z9 = 0.5*I*g2*T-0.5*I*g1*U+pi/2.0*g2-pi/2.0*g1+pi/2.0*g3;
	 Dz(0,0) = z1;
	 Dz(1,1) = z3; Dz(1,2) = -z4; Dz(1,3) = z4; Dz(1,4) = -z5;
	 Dz(2,1) = z4; Dz(2,2) = z3; Dz(2,3) = z5; Dz(2,4) = z4;
	 Dz(3,1) = -z4; Dz(3,2) = z5; Dz(3,3) = z3; Dz(3,4) = -z4;
	 Dz(4,1) = -z5; Dz(4,2) = -z4; Dz(4,3) = z4; Dz(4,4) = z3;
	 Dz(5,5) = z2;
	 Dz(6,6) = z6; Dz(6,7) = -z7; Dz(6,8) = z7; Dz(6,9) = -z5;
	 Dz(7,6) = z7; Dz(7,7) = z6; Dz(7,8) = z5; Dz(7,9) = z7;
	 Dz(8,6) = -z7; Dz(8,7) = z5; Dz(8,8) = z6; Dz(8,9) = -z7;
	 Dz(9,6) = -z5; Dz(9,7) = -z7; Dz(9,8) = z7; Dz(9,9) = z6;
	 Dz(10,10) = z8; Dz(10,11) = -z9;
	 Dz(11,10) = z9; Dz(11,11) = z8;
	 Dz(12,12) = z8; Dz(12,13) = z9;
	 Dz(13,12) = -z9; Dz(13,13) = z8;

	 G2 = Gamma2(U,T);
       }
     }
     break;
   case QQWG:
     {
       unsigned int numGauge = 1, numBrokenGauge = 6;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 1.0/sqrt(2);
       R0(1,0) = 1.0/sqrt(2);
       R0(2,0) = cos/2.0;
       R0(3,0) = sin/2.0;
       R0(4,0) = -cos/2.0;
       R0(5,0) = -sin/2.0;
       if (oneLoop) {
	 double g1 = g_Lu;
	 double g2 = g_Ld;
	 double gW = g_W;
	 
	 Complex w1 = -0.5*(T+U) + 0.75*I*pi;
	 Complex w2 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = w1;
	 Dw(2,2) = Dw(3,3) = Dw(4,4) = Dw(5,5) = w2;
         
	 Complex z1 = gW*g1*T - gW*g2*U - I*pi*(g1*g2+g1*gW-g2*gW);
	 Complex z2 = gW*g1*U - gW*g2*T - I*pi*(g2*g1+g1*gW-g2*gW);
	 Complex z3 = -I*pi*g1*g1;
	 Complex z4 = -I*pi*g2*g2;
	 Dz(0,0) = z1;
	 Dz(1,1) = z2;
	 Dz(2,2) = z3;
	 Dz(3,3) = z3;
	 Dz(4,4) = z4;
	 Dz(5,5) = z4;
         
	 G2(0,0) = -7.0/4.0*I*pi + (U+T);
       }
     }
     break;
    case QQBG:
      {
	unsigned int numGauge = 1, numBrokenGauge = 4;
	R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
	G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
	Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
	Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
	R0(0,0) = -sin;
	R0(1,0) = cos;
	R0(2,0) = -sin;
	R0(3,0) = cos;
	if (oneLoop) {
	  double g1 = g_Lu;
	  double g2 = g_Ld;
	  Complex w2 = 0.25*I*pi;
	  Dw(0,0) = Dw(1,1) = Dw(2,2) = Dw(3,3) = w2;
	  Complex z3 = -I*pi*g1*g1;
	  Complex z4 = -I*pi*g2*g2;
	  Dz(0,0) = z3;
	  Dz(1,1) = z3;
	  Dz(2,2) = z4;
	  Dz(3,3) = z4;
	  G2(0,0) = Gamma2Singlet()(0,0);
	}
      }
      break;
   case QQGG:
   case QtQtGG:
     {
       unsigned int numGauge = 3, numBrokenGauge = 6;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge); 
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = R0(3,0) = 1.0;
       R0(1,1) = R0(4,1) = 1.0;
       R0(2,2) = R0(5,2) = 1.0;
       if (oneLoop) {
	 double g1 = g_Lu;
	 double g2 = g_Ld;
	 Complex w2 = 0.25*I*pi;
	 Dw(0,0) = Dw(1,1) = Dw(2,2) = Dw(3,3) = Dw(4,4) = Dw(5,5) = w2;            
	 Complex z3 = -I*pi*g1*g1;
	 Complex z4 = -I*pi*g2*g2;
	 Dz(0,0) = Dz(1,1) = Dz(2,2) = z3;
	 Dz(3,3) = Dz(4,4) = Dz(5,5) = z4;
	 G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
	 G2 *= 0.0;
	 G2(0,0) = G2(1,1) = G2(2,2) = Gamma2Singlet()(0,0);
       }
     }
     break;
   case UUBB:
   case DDBB:
   case EEBB:
     {
       unsigned int numGauge = 1, numBrokenGauge = 4;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = sin2;
       R0(1,0) = -sin*cos;
       R0(2,0) = -sin*cos;
       R0(3,0) = cos2;
       if (oneLoop) {
	 double g1(0.);
	 if (process==UUBB) {
	   g1 = g_Ru;
	 }
	 else if (process==DDBB) {
	   g1 = g_Rd;
	 }
	 else if (process==EEBB) {
	   g1 = g_Re;
	 }
	 else
	   assert(false);
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = -I*pi*g1*g1;
	 Dz(0,0) = Dz(1,1) = Dz(2,2) = Dz(3,3) = z1;
       }
     }
     break;
   case UUPhiPhi:
   case DDPhiPhi:
   case EEPhiPhi:
     {
       unsigned int numGauge = 1, numBrokenGauge = 5;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = 1.0;
        R0(1,0) = 0.5;
        R0(2,0) = -0.5*I;
        R0(3,0) = 0.5*I;
        R0(4,0) = 0.5;
        if (oneLoop) {
	  double g1(0.);
            if (process==UUPhiPhi) {
               g1 = g_Ru;
            }
            else if (process==DDPhiPhi) {
               g1 = g_Rd;
            }
            else if (process==EEPhiPhi) {
               g1 = g_Re;
            }
	    double g3 = g_phiPlus;
            Dw(0,0) = Dw(1,4) = Dw(4,1) = 0.25*I*pi;
            Dw(2,3) = Dw(3,2) = -0.25*I*pi;
	    Complex z1 = 2.0*g3*g1*(T-U) - I*pi*(g3*g3+g1*g1);
	    Complex z2 = 0.5*I*g1*g1;
	    Complex z3 = -I*pi*g1*g1;
	    Complex z4 = 0.25*I*pi;
            Dz(0,0) = z1;
            Dz(1,1) = z3;		Dz(1,2) = -z2;		Dz(1,3) = z2;		Dz(1,4) = -z4;
            Dz(2,1) = z2;		Dz(2,2) = z3;		Dz(2,3) = z4;		Dz(2,4) = z2;
            Dz(3,1) = -z2;		Dz(3,2) = z4;		Dz(3,3) = z3;		Dz(3,4) = -z2;
            Dz(4,1) = -z4;		Dz(4,2) = -z2;		Dz(4,3) = z2;		Dz(4,4) = z3;
            G2(0,0) = Gamma2Singlet()(0,0);
	}
     }
     break;
   case UUBG:
   case DDBG:
     {
       unsigned int numGauge = 1, numBrokenGauge = 2;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       R0(0,0) = -sin;
       R0(1,0) = cos;
       if (oneLoop) {
	 double g1(0.);
	   if (process==UUBG) {
	     g1 = g_Ru;
	   }
	   else if (process==DDBG) {
	     g1 = g_Rd;
	   }
	   else
	     assert(false);
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = -I*pi*g1*g1;
	 Dz(0,0) = Dz(1,1) = z1;
       }
     }
     break;
   case UUGG:
   case tRtRGG:
   case DDGG:
     {
       unsigned int numGauge = 3, numBrokenGauge = 3;
       R0 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numGauge);
       G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
       Dw = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       Dz = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
       if (oneLoop) {
	 double g1(0.);
	 if ((process==UUGG)||(process==tRtRGG)) {
	   g1 = g_Ru;
	 }
	 else if (process==DDGG) {
	   g1 = g_Rd;
	 }
	 else
	   assert(false);
	 // There is no Dw contribution for two SU(2) singlets.
	 Complex z1 = -I*pi*g1*g1;
	 Dz(0,0) = Dz(1,1) = Dz(2,2) = z1;
       }
     }
     break;
   default:
     assert(false);
   }
   
   double aW = ElectroWeakReweighter::coupling()->aW(mu);
   double aZ = ElectroWeakReweighter::coupling()->aZ(mu);
   Energy mZ = ElectroWeakReweighter::coupling()->mZ();
   Energy mW = ElectroWeakReweighter::coupling()->mW();
   
   if (!oneLoop) {
      return R0;
   } 
   boost::numeric::ublas::matrix<Complex> temp(R0.size1(),R0.size2());
   boost::numeric::ublas::axpy_prod(R0,G2,temp);
   R0+=aW/(4.0*pi)*4.0*log(mW/mu)*temp;
   boost::numeric::ublas::axpy_prod(Dw,R0,temp);
   R0+=aW/(4.0*pi)*4.0*log(mW/mu)*temp;
   boost::numeric::ublas::axpy_prod(Dz,R0,temp);
   R0+=aZ/(4.0*pi)*4.0*log(mZ/mu)*temp;
   return R0;
}
