// -*- C++ -*-
//
// HighEnergyMatching.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#include "HighEnergyMatching.h"
#include "ElectroWeakReweighter.h"
#include "GroupInvariants.h"
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace Herwig;
using namespace HighEnergyMatching;
using namespace GroupInvariants;
using namespace EWProcess;

namespace {

/**
 *   \f$M_N\f$, this matrix is used for identical particle scattering
 */
boost::numeric::ublas::matrix<Complex> M_N(unsigned int suN) {
  double N(suN);
  boost::numeric::ublas::matrix<Complex> M(2,2);
  M(0,0) = -1.0/N;
  M(0,1) = 2.0;
  M(1,0) = 0.5 - 1.0/(2.0*N*N);
  M(1,1) = 1.0/N;
  return M;
}

#ifdef ThePEG_HAS_UNITS_CHECKING
void axpy_prod_local(const boost::numeric::ublas::matrix<Complex>  & A,
		     const boost::numeric::ublas::matrix<complex<InvEnergy2> > & B,
		     boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}

void axpy_prod_local(const boost::numeric::ublas::matrix<complex<InvEnergy2> > & A,
		       const boost::numeric::ublas::matrix<Complex> & B,
		       boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}
#else
void axpy_prod_local(const boost::numeric::ublas::matrix<Complex> & A,
		     const boost::numeric::ublas::matrix<Complex> & B,
		     boost::numeric::ublas::matrix<Complex> & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  axpy_prod(A,B,C);
}
#endif

}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
HighEnergyMatching::highEnergyMatching(Energy highScale, 
				       Energy2 s, Energy2 t, Energy2 u,
				       EWProcess::Process process,
				       bool oneLoop, bool includeAlphaS2) {
  switch (process) {
  case QQQQ:     case QQQQiden:
  case QtQtQQ:   case QQUU:
  case QtQtUU:   case QQtRtR:
  case QQDD:     case QtQtDD:
  case QQLL:     case QQEE:
  case UUUU:     case UUUUiden:
  case tRtRUU:   case UUDD:
  case tRtRDD:   case UULL:
  case UUEE:     case DDDD:
  case DDDDiden: case DDLL:
  case DDEE:     case LLLL:
  case LLLLiden: case LLEE:
  case EEEE:     case EEEEiden:
    return SpinHalfMatching(highScale,s,t,u,process,oneLoop,includeAlphaS2);
    break;         
  case QQWW:
  case QQWG:
  case QQBG:
  case QQGG:
  case QtQtGG:
  case LLWW:
  case UUBB:
  case UUBG:
  case UUGG:
  case tRtRGG:
  case DDBB:
  case DDBG:
  case DDGG:
  case EEBB:
    return Spin1HighMatching(highScale,s,t,u,process,oneLoop,includeAlphaS2);
    break;
  case QQPhiPhi:
  case LLPhiPhi:
  case UUPhiPhi:   
  case DDPhiPhi:
  case EEPhiPhi:
    return Spin0HighMatching(highScale,s,t,u,process,oneLoop,includeAlphaS2);
    break;
  default:
    assert(false);
    break;
  }
}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
HighEnergyMatching::SpinHalfMatching(Energy highScale, 
				     Energy2 s, Energy2 t, Energy2 u,
				     EWProcess::Process process, 
				     bool oneLoop, bool includeAlphaS2) {
  using Constants::pi;
  boost::numeric::ublas::matrix<complex<InvEnergy2> > highC;
  Energy Q = highScale;
  double a1 = ElectroWeakReweighter::coupling()->a1(Q);
  double a2 = ElectroWeakReweighter::coupling()->a2(Q);
  double a3 = ElectroWeakReweighter::coupling()->a3(Q);
  double y_t = ElectroWeakReweighter::coupling()->y_t(Q);
  unsigned int order = !oneLoop ? 0 : 1;
  double Yi(0.),Yf(0.);
   
  Complex Ls = MinusLog(-s/(Q*Q));
   
  double C_A_2 = C_A(2);
  double C_A_3 = C_A(3);
  double C_F_2 = C_F(2);
  double C_F_3 = C_F(3);
  double C_d_2 = C_d(2);
  double C_d_3 = C_d(3);
  double C_1_2 = C_1(2);
  double C_1_3 = C_1(3);
  Complex W = WFunction(Q,s);
  Complex X_2_st    = XNFunction(2,Q,s,t);
  //Complex X_2_su    = XNFunction(2,Q,s,u);
  Complex X_3_st    = XNFunction(3,Q,s,t);
  Complex X_3_su    = XNFunction(3,Q,s,u);
  Complex PI1       = PI1_function(Q,s);
  Complex fTilde_st = fTildeFunction(Q,s,t);
  Complex fTilde_su = fTildeFunction(Q,s,u);
   
  switch (process) {
     
  case QQQQ:
    // NOTE this 4x1 column vector highC is given by (C_11,C_21,C_12,C_22), 
    // where C_12 is the coeff. for the term that transforms under SU(2) but not SU(3)
    Yi = Yf = 1./6.;
    highC.resize(4,1);
    highC(0,0) = ZERO;
    highC(2,0) = 4.0*pi*a2 / s;
    highC(1,0) = 4.0*pi*a3 / s;
    highC(3,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-2.0*a2*a3*fTilde_st);
      highC(2,0) += (1.0/s)*(a2*a2*(X_2_st-(C_d_2+C_A_2)/4.0*fTilde_st) + 
			     2.0*(a1*a2*Yi*Yf+a2*a3*C_F_3)*W - 
			     2.0*a1*a2*Yi*Yf*fTilde_st);
      highC(1,0) += (1.0/s)*(2.0*(a1*a3*Yi*Yf+a2*a3*C_F_2)*W - 
			     2.0*a1*a3*Yi*Yf*fTilde_st);
      highC(3,0) += (1.0/s)*(-1.0*(a2*a2*C_1_2+a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     2.0*(a1*a3*Yi*Yf*C_F_3+a1*a2*Yi*Yf*C_F_2+a1*a1*Yi*Yi*Yf*Yf)*W);
      if (includeAlphaS2) {
	highC(1,0) += (1.0/s)*(a3*a3*(X_3_st-(C_d_3+C_A_3)/4.0*fTilde_st));
	
	highC(3,0) += (1.0/s)*(-1.0*(a3*a3*C_1_3)*fTilde_st);
      }
    }
    break;
     
  case QQQQiden:
    {    
      Process parentCase = QQQQ;
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCs_st = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCs_ts = highEnergyMatching(Q,t,s,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCt_st(4,1);         
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCs_ts_2x2(2,2);
      highCs_ts_2x2(0,0) = highCs_ts(0,0);
      highCs_ts_2x2(1,0) = highCs_ts(1,0);
      highCs_ts_2x2(0,1) = highCs_ts(2,0);
      highCs_ts_2x2(1,1) = highCs_ts(3,0);
      boost::numeric::ublas::matrix<Complex> temp(2,2);
      temp = boost::numeric::ublas::trans(M_N(3));
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCt_st_2x2(2,2),temp2(2,2);
      axpy_prod_local(highCs_ts_2x2,temp,temp2);
      axpy_prod_local(M_N(2),temp2,highCt_st_2x2);
      highCt_st(0,0) = highCt_st_2x2(0,0);
      highCt_st(1,0) = highCt_st_2x2(1,0);
      highCt_st(2,0) = highCt_st_2x2(0,1);
      highCt_st(3,0) = highCt_st_2x2(1,1);
      highC = highCs_st + highCt_st;
    }
    break;
  case QtQtQQ:
    {
      highC.resize(4,1);
      Process parentCase = QQQQ;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = 1.0/6.0; // Hypercharge of the non-3rd-gen doublet (still a quark doublet).
      if (order >= 1) {
	highC(2,0) += y_t*y_t*a2/(4.0*pi*s)*(3.0/2.0-0.5*Ls);
	highC(1,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0/2.0-0.5*Ls);
	highC(3,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(-5.0/12.0-1.0/12.0*Ls);
      }
    }
    break;
  case QQUU:
    Yi = 1./6.; Yf = 2./3.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a3 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*((a1*a3*(Yi*Yi+Yf*Yf)+a2*a3*C_F_2)*W + 
			     2.0*a1*a3*Yi*Yf*fTilde_su);
      highC(1,0) += (1.0/s)*(a1*a1*Yi*Yf*PI1 + 
			     (a1*a2*Yi*Yf*C_F_2+2.0*a1*a3*Yi*Yf*C_F_3+
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
      if (includeAlphaS2) {
      	highC(0,0) += (1.0/s)*(a3*a3*(X_3_su+(C_d_3-C_A_3)/4.0*fTilde_su));
      	highC(1,0) += (1.0/s)*((a3*a3*C_1_3+a1*a1*Yi*Yi*Yf*Yf)*fTilde_su);
      }
    }
    break;
  case QtQtUU:
    {
      highC.resize(2,1);
      Process parentCase = QQUU;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = 2./3.;
      if (order >= 1) {
	highC(0,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0/2.0-0.5*Ls);
	highC(1,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(-5.0/12.0-1.0/12.0*Ls);
      }
    }
    break;         
  case QQtRtR:
    {
      highC.resize(2,1);
      Process parentCase = QQUU;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = 1.0/6.0;
      if (order >= 1) {
	highC(0,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0-Ls);
	highC(1,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(5.0/3.0-2.0/3.0*Ls);
      }
    }
    break;
  case QQDD:
    Yi = 1./6.; Yf = -1./3.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a3 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*((a1*a3*(Yi*Yi+Yf*Yf)+a2*a3*C_F_2)*W + 
			     2.0*a1*a3*Yi*Yf*fTilde_su);
      highC(1,0) += (1.0/s)*(a1*a1*Yi*Yf*PI1 + 
			     (a1*a2*Yi*Yf*C_F_2+2.0*a1*a3*Yi*Yf*C_F_3+
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
      if (includeAlphaS2) {
      	highC(0,0) += (1.0/s)*(a3*a3*(X_3_su+(C_d_3-C_A_3)/4.0*fTilde_su));
      	highC(1,0) += (1.0/s)*((a3*a3*C_1_3+a1*a1*Yi*Yi*Yf*Yf)*fTilde_su);
      }
    }
    break;
  case QtQtDD:
    {
      highC.resize(2,1);
      Process parentCase = QQDD;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = -1./3.;
      if (order >= 1) {
	highC(0,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0/2.0-0.5*Ls);
	highC(1,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(-5.0/12.0-1.0/12.0*Ls);
      }
    }
    break;
  case QQLL:
    Yi = 1./6.; Yf = -1./2.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a2 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a2*a2*(X_2_st-(C_d_2+C_A_2)/4.0*fTilde_st) + 
			     (a2*a3*C_F_3 + a1*a2*(Yi*Yi+Yf*Yf))*W - 
			     2.0*a1*a2*Yi*Yf*fTilde_st);
      highC(1,0) += (1.0/s)*(-1.0*(a2*a2*C_1_2+a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3+2.0*a1*a2*Yi*Yf*C_F_2+
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;
  case QQEE:
    Yi = 1./6.; Yf = -1.;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a1*a1*Yi*Yi*Yf*Yf*fTilde_su + a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3 + a1*a2*Yi*Yf*C_F_2 + 
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;
  case UUUU:
    Yi = Yf = 2./3.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a3 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-2.0*a1*a3*Yi*Yf*fTilde_st +
			     a1*a3*(Yi*Yi+Yf*Yf)*W);
      highC(1,0) += (1.0/s)*(-1.0*(a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     (2.0*a1*a3*Yi*Yf*C_F_3+a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
      if (includeAlphaS2) {
         highC(0,0) += (1.0/s)*(a3*a3*(X_3_st-(C_d_3+C_A_3)/4.0*fTilde_st));
         highC(1,0) += (1.0/s)*(-1.0*(a3*a3*C_1_3)*fTilde_st);
      }
    }
    break;
  case UUUUiden:
    {
      Process parentCase = UUUU;
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCs_st = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> >
	highCs_ts = highEnergyMatching(Q,t,s,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCt_st;
	axpy_prod_local(M_N(3),highCs_ts,highCt_st);
      highC = highCs_st + highCt_st;
    }
    break;
  case tRtRUU:
    {
      highC.resize(2,1);
      Process parentCase = UUUU;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = 2./3.;
      if (order >= 1) {
	highC(0,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0-Ls);
	highC(1,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(5.0/3.0-2.0/3.0*Ls);
      }
    }
    break;
  case UUDD:
    Yi = 2./3.; Yf = -1./3.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a3 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-2.0*a1*a3*Yi*Yf*fTilde_st +
			     a1*a3*(Yi*Yi+Yf*Yf)*W);
      highC(1,0) += (1.0/s)*(-1.0*(a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     (2.0*a1*a3*Yi*Yf*C_F_3+a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
      if (includeAlphaS2) {
	highC(0,0) += (1.0/s)*(a3*a3*(X_3_st-(C_d_3+C_A_3)/4.0*fTilde_st));
	highC(1,0) += (1.0/s)*(-1.0*(a3*a3*C_1_3)*fTilde_st);
      }
    }
    break;
  case tRtRDD:
    {
      highC.resize(2,1);
      Process parentCase = UUDD;
      highC = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      double Y = -1./3.;
      if (order >= 1) {
	highC(0,0) += y_t*y_t*a3/(4.0*pi*s)*(1.0-Ls);
	highC(1,0) += y_t*y_t*a1*Y/(4.0*pi*s)*(5.0/3.0-2.0/3.0*Ls);
      }
    }
    break;
  case UULL:
    Yi = 2./3.; Yf = -0.5;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a1*a1*Yi*Yi*Yf*Yf*fTilde_su + a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3 + a1*a2*Yi*Yf*C_F_2 + 
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;
  case UUEE:
    Yi = 2./3.; Yf = -1.;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-1.0*a1*a1*Yi*Yi*Yf*Yf*fTilde_st + a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3 + 
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;
  case DDDD:
    Yi = Yf = -1./3.;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a3 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-2.0*a1*a3*Yi*Yf*fTilde_st +
			     a1*a3*(Yi*Yi+Yf*Yf)*W);
      highC(1,0) += (1.0/s)*(-1.0*(a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     (2.0*a1*a3*Yi*Yf*C_F_3+a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
      if (includeAlphaS2) {
         highC(0,0) += (1.0/s)*(a3*a3*(X_3_st-(C_d_3+C_A_3)/4.0*fTilde_st));
         highC(1,0) += (1.0/s)*(-1.0*(a3*a3*C_1_3)*fTilde_st);
      }
    }
    break;
  case DDDDiden:
    {
      Process parentCase = DDDD;
      boost::numeric::ublas::matrix<complex<InvEnergy2> >
	highCs_st = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> >
	highCs_ts = highEnergyMatching(Q,t,s,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCt_st;
      axpy_prod_local(M_N(3),highCs_ts,highCt_st);
      highC = highCs_st + highCt_st;
    }
    break;
  case DDLL:
    Yi = -1./3.; Yf = -0.5;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a1*a1*Yi*Yi*Yf*Yf*fTilde_su + a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3 + a1*a2*Yi*Yf*C_F_2 + 
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;         
  case DDEE:
    Yi = -1./3.; Yf = -1.;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-1.0*a1*a1*Yi*Yi*Yf*Yf*fTilde_st + a1*a1*Yi*Yf*PI1 + 
			     (a1*a3*Yi*Yf*C_F_3 + 
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;         
  case LLLL:
    Yi = Yf = -0.5;
    highC.resize(2,1);
    highC(0,0) = 4.0*pi*a2 / s;
    highC(1,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a2*a2*(X_2_st-(C_d_2+C_A_2)/4.0*fTilde_st) + 
			     2.0*a1*a2*Yi*Yf*W -
			     2.0*a1*a2*Yi*Yf*fTilde_st);
      highC(1,0) += (1.0/s)*(-1.0*(a2*a2*C_1_2+a1*a1*Yi*Yi*Yf*Yf)*fTilde_st + 
			     a1*a1*Yi*Yf*PI1 + 
			     2.0*(a1*a2*Yi*Yf*C_F_2+a1*a1*Yi*Yi*Yf*Yf)*W);
    }
    break;
  case LLLLiden:
    {
      Process parentCase = LLLL;
      boost::numeric::ublas::matrix<complex<InvEnergy2> >
  	highCs_st = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> >
  	highCs_ts = highEnergyMatching(Q,t,s,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > highCt_st;
      axpy_prod_local(M_N(2), highCs_ts, highCt_st);
      highC = highCs_st + highCt_st;
    }
    break;
  case LLEE:
    Yi = -0.5; Yf = -1.;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(a1*a1*Yi*Yi*Yf*Yf*fTilde_su + a1*a1*Yi*Yf*PI1 + 
			     (a1*a2*Yi*Yf*C_F_2 +  
			      a1*a1*(Yi*Yi*Yi*Yf+Yf*Yf*Yf*Yi))*W);
    }
    break;
  case EEEE:
    Yi = Yf = -1.;
    highC.resize(1,1);
    highC(0,0) = 4.0*pi*a1*Yi*Yf / s;
    if (order >= 1) {
      highC(0,0) += (1.0/s)*(-1.0*a1*a1*Yi*Yi*Yf*Yf*fTilde_st + a1*a1*Yi*Yf*PI1 + 
			     2.0*a1*a1*Yi*Yi*Yf*Yf*W);
    }
    break;
  case EEEEiden: 
    {
      Process parentCase = EEEE;
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCs_st = highEnergyMatching(Q,s,t,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCs_ts = highEnergyMatching(Q,t,s,u,parentCase,oneLoop,includeAlphaS2);
      boost::numeric::ublas::matrix<complex<InvEnergy2> > 
	highCt_st = highCs_ts;
      highC = highCs_st + highCt_st;
    }
    break;
  default:
    assert(false);
  }
  return highC;
}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
HighEnergyMatching::Spin1HighMatching(Energy highScale, 
				      Energy2 s, Energy2 t, Energy2 u,
				      EWProcess::Process process,
				      bool oneLoop, bool includeAlphaS2) {
  using Constants::pi;
   
  unsigned int order = !oneLoop ? 0 : 1;
  // (If crossed graphs, swap s and t here)
  Complex L_s = MinusLog(-s/(highScale*highScale));
  Complex L_t = MinusLog(-t/(highScale*highScale));
  Complex L_u = MinusLog(-u/(highScale*highScale));
  Complex L_s2 = L_s*L_s;
  Complex L_t2 = L_t*L_t;
  Complex L_u2 = L_u*L_u;
   
  // Tree-Level:
  // Topology types defined. Note each of these is a vector of 5 entries. They are the coefficients 
  // for the dirac structures M_0, M_1, M_4, M_5, and M_6 for vector boson production.
  boost::numeric::ublas::vector<complex<InvEnergy2> > R1(5);
  for(unsigned int ix=0;ix<5;++ix) R1[ix] = ZERO;
  R1[0] = -1.0/t;
  R1[1] = -2.0/t;
  R1[2] = ZERO;
  R1[3] = ZERO;
  R1[4] = ZERO;
  boost::numeric::ublas::vector<complex<InvEnergy2> > R1_bar(5);
  for(unsigned int ix=0;ix<5;++ix) R1_bar[ix] = ZERO;
  R1_bar[0] = -1.0/u;
  boost::numeric::ublas::vector<complex<InvEnergy2> > R2(5);
  for(unsigned int ix=0;ix<5;++ix) R2[ix] = ZERO;
  R2[1] = -1.0/s*2.0;
  // Topologies T1:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1a(5);
  for(unsigned int ix=0;ix<5;++ix) T1a[ix] = ZERO;
  T1a[0] = 1.0/(t*u)*(-3.0*t*L_s2-(s+4.0*t)*L_t2+2.0*(s+4.0*t)*L_s*L_t+2.0*u*L_t-
                      pi*pi*(7.0/6.0*s+25.0/6.0*t)-4.0*u);
  T1a[1] = 1.0/(u*u*t*s)*(0.5*t*(9.0*s*s+14.0*s*t+7.0*t*t)*L_s2+s*(2.0*s+t)*(s+2.0*t)*L_t2-
                          2.0*(2.0*s*s*s+9.0*s*s*t+10.0*s*t*t+4.0*t*t*t)*L_s*L_t-
                          2.0*t*t*u*L_s-2.0*u*s*(2.0*s+3.0*t)*L_t+
                          pi*pi*(7.0/3.0*s*s*s+125.0/12.0*s*s*t+71.0/6.0*s*t*t+
                                 19.0/4.0*t*t*t)-
                          8.0*s*s*s-20.0*s*s*t-16.0*s*t*t-4.0*t*t*t);
  T1a[2] = 1.0/(t*u*u)*(-t*(3.0*s+4.0*t)*L_s2-(s*s+5.0*s*t+5.0*t*t)*L_t2+
                        2.0*t*(3.0*s+4.0*t)*L_s*L_t+2.0*u*t*(2.0*s+t)*L_s/s-
                        2.0*u*t*L_t+pi*pi*(s*s/6.0-8.0/3.0*s*t-23.0/6.0*t*t)+
                        4.0*t*t*t/s+4.0*s*t+8.0*t*t);
  T1a[3] = T1a[2];
  T1a[4] = GeV2/(t*u*u*u)*(-4.0*t*(s+2.0*t)*(L_s-L_t)*(L_s-L_t)+
			   4.0*u*(3.0*s+5.0*t)*(L_s-L_t)-4.0*pi*pi*t*(s+2.0*t)-4.0*u*u);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1b(5);
  for(unsigned int ix=0;ix<5;++ix) T1b[ix] = ZERO;
  T1b[0] = 1.0/(t*u*s*s)*(-s*t*(2.0*s+3.0*t)*L_u2+s*u*(s+3.0*t)*L_t2+
			  2.0*s*(s*s+3.0*s*t+3.0*t*t)*L_u*L_t+s*s*t*L_u+s*s*u*L_t-
			  pi*pi*(7.0/6.0*s*s*s+3.0*s*s*t+3.0*s*t*t)+2.0*s*s*s);
  T1b[1] = 1.0/(t*s*s*u)*(3.0*s*t*u*L_u2+s*u*(2.0*s+3.0*t)*L_t2-
			  2.0*s*u*(2.0*s+3.0*t)*L_u*L_t+2.0*s*s*u*L_t-
			  pi*pi*(7.0/3.0*s*s*s+16.0/3.0*s*s*t+3.0*s*t*t)+
			  4.0*s*s*s+4.0*s*s*t);
  T1b[2] = 1.0/(t*u*s*s)*(-3.0*s*t*u*(L_u-L_t)*(L_u-L_t)+4.0*s*s*t*L_u+4.0*s*s*u*L_t+
			  pi*pi*(3.0*s*s*t+3.0*s*t*t)+8.0*s*s*s);
  T1b[3] = 1.0/(t*u*s*s)*(s*t*(2.0*s+3.0*t)*L_u2-s*u*(s+3.0*t)*L_t2+6.0*s*t*u*L_u*L_t+
			  pi*pi*(-1.0/6.0*s*s*s+3.0*s*s*t+3.0*s*t*t));
  T1b[4] = 12.0*GeV2/(t*u)*(L_t-L_u);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1c(5);
  for(unsigned int ix=0;ix<5;++ix) T1c[ix] = ZERO;
  T1c[0] = 1.0/t*(2.0*L_s*L_t-7.0*pi*pi/6.0-L_t2);
  T1c[1] = 1.0/(t*u*u)*(s*t*L_s2-(2.0*s*s+3.0*s*t+2.0*t*t)*L_t2+
			2.0*(2.0*s*s+3.0*s*t+2.0*t*t)*L_s*L_t+2.0*t*u*(L_s-L_t)-
			pi*pi*(7.0/3.0*s*s+11.0/3.0*s*t+7.0/3.0*t*t));
  T1c[2] = 1.0/(t*u*u)*(t*(3.0*s+2.0*t)*(L_s-L_t)*(L_s-L_t)+2.0*u*t*L_s+
			2.0*u*(2.0*s+t)*L_t+pi*pi*t*(3.0*s+2.0*t)+8.0*u*u);
  T1c[3] = T1c[2];
  T1c[4] = GeV2/(t*u*u*u)*(4.0*t*(2.0*s+t)*(L_s-L_t)*(L_s-L_t)-4.0*u*(3.0*s+t)*(L_s-L_t)+
			   4.0*pi*pi*t*(2.0*s+t)-4.0*u*u);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1d(5);
  for(unsigned int ix=0;ix<5;++ix) T1d[ix] = ZERO;
  T1d[2] = 1.0/s*(-2.0*L_s+4.0);
  T1d[3] = T1d[2];
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1a_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T1a_bar[ix] = ZERO;
  T1a_bar[0] = 1.0/(u*t)*(-3.0*u*L_s2-(s+4.0*u)*L_u2+2.0*(s+4.0*u)*L_s*L_u+2.0*t*L_u-
			  pi*pi*(7.0/6.0*s+25.0/6.0*u)-4.0*t);
  T1a_bar[1] = 2.0*T1a_bar[0] - 
    1.0/(t*t*u*s)*(0.5*u*(9.0*s*s+14.0*s*u+7.0*u*u)*L_s2+s*(2.0*s+u)*(s+2.0*u)*L_u2-
		   2.0*(2.0*s*s*s+9.0*s*s*u+10.0*s*u*u+4.0*u*u*u)*L_s*L_u-
		   2.0*u*u*t*L_s-2.0*t*s*(2.0*s+3.0*u)*L_u+
		   pi*pi*(7.0/3.0*s*s*s+125.0/12.0*s*s*u+71.0/6.0*s*u*u+
			  19.0/4.0*u*u*u)-
		   8.0*s*s*s-20.0*s*s*u-16.0*s*u*u-4.0*u*u*u);
  T1a_bar[2] = 1.0/(u*t*t)*(-u*(3.0*s+4.0*u)*L_s2-(s*s+5.0*s*u+5.0*u*u)*L_u2+
			    2.0*u*(3.0*s+4.0*u)*L_s*L_u+2.0*t*u*(2.0*s+u)*L_s/s-
			    2.0*t*u*L_u+pi*pi*(s*s/6.0-8.0/3.0*s*u-23.0/6.0*u*u)+
			    4.0*u*u*u/s+4.0*s*u+8.0*u*u);
  T1a_bar[3] = T1a_bar[2];
  T1a_bar[4] = -GeV2/(u*t*t*t)*(-4.0*u*(s+2.0*u)*(L_s-L_u)*(L_s-L_u)+
				4.0*t*(3.0*s+5.0*u)*(L_s-L_u)-4.0*pi*pi*u*(s+2.0*u)-4.0*t*t);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1b_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T1b_bar[ix] = ZERO;
  T1b_bar[0] = 1.0/(u*t*s*s)*(-s*u*(2.0*s+3.0*u)*L_t2+s*t*(s+3.0*u)*L_u2+
			      2.0*s*(s*s+3.0*s*u+3.0*u*u)*L_t*L_u+
			      s*s*u*L_t+s*s*t*L_u-
			      pi*pi*(7.0/6.0*s*s*s+3.0*s*s*u+3.0*s*u*u)+2.0*s*s*s);
  T1b_bar[1] = 2.0*T1b_bar[0] - 
    1.0/(u*s*s*t)*(3.0*s*u*t*L_t2+s*t*(2.0*s+3.0*u)*L_u2-
		   2.0*s*t*(2.0*s+3.0*u)*L_t*L_u+2.0*s*s*t*L_u-
		   pi*pi*(7.0/3.0*s*s*s+16.0/3.0*s*s*u+3.0*s*u*u)+
		   4.0*s*s*s+4.0*s*s*u);
  T1b_bar[3] = 1.0/(u*t*s*s)*(-3.0*s*u*t*(L_t-L_u)*(L_t-L_u)+4.0*s*s*u*L_t+4.0*s*s*t*L_u+
			      pi*pi*(3.0*s*s*u+3.0*s*u*u)+8.0*s*s*s);
  T1b_bar[2] = 1.0/(u*t*s*s)*(s*u*(2.0*s+3.0*u)*L_t2-s*t*(s+3.0*u)*L_u2+6.0*s*u*t*L_t*L_u+
			      pi*pi*(-1.0/6.0*s*s*s+3.0*s*s*u+3.0*s*u*u));
  T1b_bar[4] = -12.0*GeV2/(u*t)*(L_u-L_t);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T1c_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T1c_bar[ix] = ZERO;
  T1c_bar[0] = 1.0/u*(2.0*L_s*L_u-7.0*pi*pi/6.0-L_u2);
  T1c_bar[1] = 2.0*T1c_bar[0] - 
    1.0/(u*t*t)*(s*u*L_s2-(2.0*s*s+3.0*s*u+2.0*u*u)*L_u2+
		 2.0*(2.0*s*s+3.0*s*u+2.0*u*u)*L_s*L_u+2.0*u*t*(L_s-L_u)-
		 pi*pi*(7.0/3.0*s*s+11.0/3.0*s*u+7.0/3.0*u*u));
  T1c_bar[2] = 1.0/(u*t*t)*(u*(3.0*s+2.0*u)*(L_s-L_u)*(L_s-L_u)+2.0*t*u*L_s+
                            2.0*t*(2.0*s+u)*L_u+pi*pi*u*(3.0*s+2.0*u)+8.0*t*t);
  T1c_bar[3] = T1c_bar[2];
  T1c_bar[4] = -GeV2/(u*t*t*t)*(4.0*u*(2.0*s+u)*(L_s-L_u)*(L_s-L_u)-4.0*t*(3.0*s+u)*(L_s-L_u)+
                               4.0*pi*pi*u*(2.0*s+u)-4.0*t*t);
  // Topologies T2:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T2a(5);
  for(unsigned int ix=0;ix<5;++ix) T2a[ix] = ZERO;
  T2a[1] = 1.0/s*(L_s/6.0-11.0/18.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T2b(5);
  for(unsigned int ix=0;ix<5;++ix) T2b[ix] = ZERO;
  T2b[1] = 1.0/s*(-2.0/3.0*L_s+22.0/9.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T2c(5);
  for(unsigned int ix=0;ix<5;++ix) T2c[ix] = ZERO;
  T2c[1] = 1.0/s*(3.0/2.0*L_s2-17.0/2.0*L_s-pi*pi/4.0+95.0/6.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T2d(5);
  for(unsigned int ix=0;ix<5;++ix) T2d[ix] = ZERO;
  T2d[1] = 1.0/s*(-4.0/3.0*L_s+14.0/9.0);
  // Topologies T3:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T3a(5);
  for(unsigned int ix=0;ix<5;++ix) T3a[ix] = ZERO;
  T3a[0] = 1.0/t*(L_t2-pi*pi/6.0);
  T3a[1] = 2.0*T3a[0];
  T3a[3] = 1.0/t*(-L_t2+pi*pi/6.0+2.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T3b(5);
  for(unsigned int ix=0;ix<5;++ix) T3b[ix] = ZERO;
  T3b[0] = 1.0/t*(-L_t+4.0);
  T3b[1] = 2.0*T3b[0];
  T3b[3] = 1.0/t*(4.0*L_t-10.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T3a_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T3a_bar[ix] = ZERO;
  T3a_bar[0] = 1.0/u*(L_u2-pi*pi/6.0);
  T3a_bar[2] = 1.0/u*(-L_u2+pi*pi/6.0+2.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T3b_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T3b_bar[ix] = ZERO;
  T3b_bar[0] = 1.0/u*(-L_u+4.0);
  T3b_bar[2] = 1.0/u*(4.0*L_u-10.0);
  // Topologies T4:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T4a(5);
  for(unsigned int ix=0;ix<5;++ix) T4a[ix] = ZERO;
  T4a[0] = 1.0/t*(L_t2-pi*pi/6.0);
  T4a[1] = 2.0*T4a[0];
  T4a[2] = 1.0/t*(-L_t2+pi*pi/6.0+2.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T4b(5);
  for(unsigned int ix=0;ix<5;++ix) T4b[ix] = ZERO;
  T4b[0] = 1.0/t*(-L_t+4.0);
  T4b[1] = 2.0*T4b[0];
  T4b[2] = 1.0/t*(4.0*L_t-10.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T4a_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T4a_bar[ix] = ZERO;
  T4a_bar[0] = 1.0/u*(L_u2-pi*pi/6.0);
  T4a_bar[3] = 1.0/u*(-L_u2+pi*pi/6.0+2.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T4b_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T4b_bar[ix] = ZERO;
  T4b_bar[0] = 1.0/u*(-L_u+4.0);
  T4b_bar[3] = 1.0/u*(4.0*L_u-10.0);
  // Topologies T5:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T5a(5);
  for(unsigned int ix=0;ix<5;++ix) T5a[ix] = ZERO;
  T5a[1] = 1.0/s*(2.0*L_s-4.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T5b(5);
  for(unsigned int ix=0;ix<5;++ix) T5b[ix] = ZERO;
  T5b[1] = 1.0/s*(-2.0*L_s2+6.0*L_s-16.0+pi*pi/3.0);
  // Topologies T6:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T6a(5);
  for(unsigned int ix=0;ix<5;++ix) T6a[ix] = ZERO;
  T6a[1] = 1.0/s*(-19.0/6.0*L_s+58.0/9.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T6b(5);
  for(unsigned int ix=0;ix<5;++ix) T6b[ix] = ZERO;
  T6b[1] = 1.0/s*(-1.0/6.0*L_s+4.0/9.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T6c(5);
  for(unsigned int ix=0;ix<5;++ix) T6c[ix] = ZERO;
  T6c[1] = 1.0/s*(2.0/3.0*L_s-16.0/9.0);
  boost::numeric::ublas::vector<complex<InvEnergy2> > T6d(5);
  for(unsigned int ix=0;ix<5;++ix) T6d[ix] = ZERO;
  T6d[1] = 1.0/s*(4.0/3.0*L_s-20.0/9.0);
  // Topology T7:
  boost::numeric::ublas::vector<complex<InvEnergy2> > T7(5);
  for(unsigned int ix=0;ix<5;++ix) T7[ix] = ZERO;
  T7[0] = 1.0/t*(-L_t+1.0);
  T7[1] = 2.0*T7[0];
  boost::numeric::ublas::vector<complex<InvEnergy2> > T7_bar(5);
  for(unsigned int ix=0;ix<5;++ix) T7_bar[ix] = ZERO;
  T7_bar[0] = 1.0/u*(-L_u+1.0);
  // Group Theory Factors / SM parameters needed for matrix elements:
  double a1 = ElectroWeakReweighter::coupling()->a1(highScale);
  double a2 = ElectroWeakReweighter::coupling()->a2(highScale);
  double a3 = ElectroWeakReweighter::coupling()->a3(highScale);
  double y_t = ElectroWeakReweighter::coupling()->y_t(highScale);
  // Traces over complex scalars and weyl fermions.
  double T_CS_3 = 0.0;
  double T_CS_2 = 0.5;
  //double T_CS_1 = 0.5;
  double T_WF_3 = 2.0*3.0;
  double T_WF_2 = 2.0*3.0;
  //double T_WF_1 = 10.0/3.0*3.0;
  double C_A_3 = 3.0;
  double C_A_2 = 2.0;
  double C_A_1 = 0.0;
  double C_F_3 = 4.0/3.0;
  double C_F_2 = 3.0/4.0;
  double C_F_1 = 1.0;
  // This is the coefficient of the delta term in G_TT
  double G_TT = 0.5;
  // This is the coeffidient of d^ABC in G_TT (non-zero for SU(3))
  double G_TT_3_D = 0.25*C_A_3; 
  double G_f = 1.0;
  // Factors TBD after fermion helicity is specified:
  double Y_Q(0.), G_Plus_U1(0.);
  double G_Plus_SU2 = 0.25;
  double G_Plus_SU3 = 1./6.;
  double G_Plus_SU3_D = 0.5;
  double Lambda_Q(0.);
  // the _s and _ew are the alpha3 and alpha1/2 parts of Lambda_Q
  double Lambda_Q_s(0.); 
  double Lambda_Q_ew(0.);
  double rho12(0.), rho13(0.);
  double rho23 = sqrt(a2*a3);
  double tRorQ = 1.0;
  boost::numeric::ublas::matrix<complex<InvEnergy2> > highC(1,1);
  switch (process) {
  case QQWW: case LLWW:
    {
      // Finish Group Theory Factors:
      if (process==QQWW) {
	Y_Q = 1.0/6.0;
	G_Plus_U1 = Y_Q*Y_Q;
	Lambda_Q = C_F_3*a3 + C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
	rho12 = Y_Q*sqrt(a1*a2);
      }
      else if (process==LLWW) {
	Y_Q = -1.0/2.0;
	G_Plus_U1 = Y_Q*Y_Q;
	Lambda_Q = C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
	rho12 = Y_Q*sqrt(a1*a2);
      }
      highC.resize(5,5);
      for (int i=0; i<5; i++) {
	highC(0,i) = G_Plus_SU2*(4.0*pi*a2)*(R1[i]+R1_bar[i]);
	highC(1,i) = G_f*4.0*pi*a2*(-0.5*R1[i]+0.5*R1_bar[i]-R2[i]);
	highC(2,i) = rho12*4.0*pi*(R1[i]+R1_bar[i]);
	highC(3,i) = rho12*4.0*pi*(R1[i]+R1_bar[i]);
	highC(4,i) = G_Plus_U1*(4.0*pi*a1)*(R1[i]+R1_bar[i]);
	if (order>=1) {
	  highC(0,i) += G_Plus_SU2*((-0.5*a2*a2*C_A_2)*(T1b[i]+T1b_bar[i])+
				    a2*(Lambda_Q-a2*C_A_2)*(T1c[i]+T1c_bar[i])+
				    0.5*a2*a2*C_A_2*(T3a[i]+T3a_bar[i])+
				    a2*(Lambda_Q-0.5*a2*C_A_2)*(T3b[i]+T3b_bar[i])+
				    0.5*a2*a2*C_A_2*(T4a[i]+T4a_bar[i])+
				    a2*(Lambda_Q-0.5*a2*C_A_2)*(T4b[i]+T4b_bar[i])+
				    a2*Lambda_Q*(T7[i]+T7_bar[i])) + 
	    G_TT*(-a2*a2*(T1a[i]+T1a_bar[i])+a2*a2*(T1b[i]+T1b_bar[i])+
		  a2*a2*(T1c[i]+T1c_bar[i])+2.0*a2*a2*T1d[i]);
	  highC(1,i) += G_f*(0.25*a2*a2*C_A_2*(T1a[i]-T1a_bar[i])+
			     a2*(0.25*a2*C_A_2-0.5*Lambda_Q)*(T1c[i]-T1c_bar[i])+
			     0.5*a2*a2*C_A_2*T2a[i]+a2*a2*T_CS_2*T2b[i]-
			     0.5*a2*a2*C_A_2*T2c[i]+a2*a2*T_WF_2*T2d[i]-
			     0.25*a2*a2*C_A_2*(T3a[i]-T3a_bar[i])-
			     0.5*a2*(Lambda_Q-0.5*a2*C_A_2)*(T3b[i]-T3b_bar[i])-
			     0.25*a2*a2*C_A_2*(T4a[i]-T4a_bar[i])-
			     0.5*a2*(Lambda_Q-0.5*a2*C_A_2)*(T4b[i]-T4b_bar[i])+
			     0.5*a2*a2*C_A_2*T5a[i]+a2*(Lambda_Q-0.5*a2*C_A_2)*T5b[i]+
			     a2*a2*C_A_2*T6a[i]+a2*a2*C_A_2*T6b[i]+
			     a2*a2*T_CS_2*T6c[i]+a2*a2*T_WF_2*T6d[i]-
			     0.5*a2*Lambda_Q*(T7[i]-T7_bar[i]));
	  highC(2,i) += rho12*(-0.5*a1*C_A_1*T1b[i]-0.5*a2*C_A_2*T1b_bar[i]+
			       (Lambda_Q-0.5*a1*C_A_1-0.5*a2*C_A_2)*(T1c[i]+T1c_bar[i])+
			       0.5*a1*C_A_1*T3a[i]+0.5*a2*C_A_2*T3a_bar[i]+
			       (Lambda_Q-0.5*a1*C_A_1)*T3b[i]+(Lambda_Q-0.5*a2*C_A_2)*T3b_bar[i]+
			       0.5*a2*C_A_2*T4a[i]+0.5*a1*C_A_1*T4a_bar[i]+
			       (Lambda_Q-0.5*a2*C_A_2)*T4b[i]+(Lambda_Q-0.5*a1*C_A_1)*T4b_bar[i]+
			       Lambda_Q*(T7[i]+T7_bar[i]));
	  highC(3,i) += rho12*(-0.5*a2*C_A_2*T1b[i]-0.5*a1*C_A_1*T1b_bar[i]+
			       (Lambda_Q-0.5*a2*C_A_2-0.5*a1*C_A_1)*(T1c[i]+T1c_bar[i])+
			       0.5*a2*C_A_2*T3a[i]+0.5*a1*C_A_1*T3a_bar[i]+
			       (Lambda_Q-0.5*a2*C_A_2)*T3b[i]+(Lambda_Q-0.5*a1*C_A_1)*T3b_bar[i]+
			       0.5*a1*C_A_1*T4a[i]+0.5*a2*C_A_2*T4a_bar[i]+
			       (Lambda_Q-0.5*a1*C_A_1)*T4b[i]+(Lambda_Q-0.5*a2*C_A_2)*T4b_bar[i]+
			       Lambda_Q*(T7[i]+T7_bar[i]));
	  highC(4,i) += G_Plus_U1*((-0.5*a1*a1*C_A_1)*(T1b[i]+T1b_bar[i])+
				   a1*(Lambda_Q-a1*C_A_1)*(T1c[i]+T1c_bar[i])+
				   0.5*a1*a1*C_A_1*(T3a[i]+T3a_bar[i])+
				   a1*(Lambda_Q-0.5*a1*C_A_1)*(T3b[i]+T3b_bar[i])+
				   0.5*a1*a1*C_A_1*(T4a[i]+T4a_bar[i])+
				   a1*(Lambda_Q-0.5*a1*C_A_1)*(T4b[i]+T4b_bar[i])+
				   a1*Lambda_Q*(T7[i]+T7_bar[i]));
	}
      }
    }
    break;
  case UUBB: case DDBB: case EEBB: 
    {
      // Finish Group Theory Factors:
      if (process==UUBB) {
	Y_Q = 2.0/3.0;
	G_Plus_U1 = Y_Q*Y_Q;
	Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
      }
      else if (process==DDBB) {
	Y_Q = -1.0/3.0;
	G_Plus_U1 = Y_Q*Y_Q;
	Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
      }
      else if (process==EEBB) {
	Y_Q = -1.0;
	G_Plus_U1 = Y_Q*Y_Q;
	Lambda_Q = Y_Q*Y_Q*C_F_1*a1;
      }
      highC.resize(1,5);
      for (int i=0; i<5; i++) {
	highC(0,i) = G_Plus_U1*(4.0*pi*a1)*(R1[i]+R1_bar[i]);
	if (order>=1) {
	  highC(0,i) += G_Plus_U1*((-0.5*a1*a1*C_A_1)*(T1b[i]+T1b_bar[i])+
				   a1*(Lambda_Q-a1*C_A_1)*(T1c[i]+T1c_bar[i])+
				   0.5*a1*a1*C_A_1*(T3a[i]+T3a_bar[i])+
				   a1*(Lambda_Q-0.5*a1*C_A_1)*(T3b[i]+T3b_bar[i])+
				   0.5*a1*a1*C_A_1*(T4a[i]+T4a_bar[i])+
				   a1*(Lambda_Q-0.5*a1*C_A_1)*(T4b[i]+T4b_bar[i])+
				   a1*Lambda_Q*(T7[i]+T7_bar[i]));
	}
      }
    }
    break;
  case QQWG:
    {
      // Finish Group Theory Factors:
      Y_Q = 1./6.;
      Lambda_Q = C_F_3*a3 + C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
      
      highC.resize(1,5);
      
      for (int i=0; i<5; i++) {
	highC(0,i) = rho23*4.0*pi*(R1[i]+R1_bar[i]);
        
	if (order>=1) {
	  
	  highC(0,i) += rho23*(-0.5*a3*C_A_3*T1b[i]-0.5*a2*C_A_2*T1b_bar[i]+
			       (Lambda_Q-0.5*a3*C_A_3-0.5*a2*C_A_2)*(T1c[i]+T1c_bar[i])+
			       0.5*a3*C_A_3*T3a[i]+0.5*a2*C_A_2*T3a_bar[i]+
			       (Lambda_Q-0.5*a3*C_A_3)*T3b[i]+(Lambda_Q-0.5*a2*C_A_2)*T3b_bar[i]+
			       0.5*a2*C_A_2*T4a[i]+0.5*a3*C_A_3*T4a_bar[i]+
			       (Lambda_Q-0.5*a2*C_A_2)*T4b[i]+(Lambda_Q-0.5*a3*C_A_3)*T4b_bar[i]+
			       Lambda_Q*(T7[i]+T7_bar[i]));
	}
      }
    }
    break;
  case QQBG:
    {
      // Finish Group Theory Factors:
      Y_Q = 1.0/6.0;
      Lambda_Q = C_F_3*a3 + C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
      rho13 = Y_Q*sqrt(a1*a3);
      
      highC.resize(1,5);
      
      for (int i=0; i<5; i++) {
	highC(0,i) = rho13*4.0*pi*(R1[i]+R1_bar[i]);
        
	if (order>=1) {
	  
	  highC(0,i) += rho13*(-0.5*a3*C_A_3*T1b[i]-0.5*a1*C_A_1*T1b_bar[i]+
			       (Lambda_Q-0.5*a3*C_A_3-0.5*a1*C_A_1)*(T1c[i]+T1c_bar[i])+
			       0.5*a3*C_A_3*T3a[i]+0.5*a1*C_A_1*T3a_bar[i]+
			       (Lambda_Q-0.5*a3*C_A_3)*T3b[i]+(Lambda_Q-0.5*a1*C_A_1)*T3b_bar[i]+
			       0.5*a1*C_A_1*T4a[i]+0.5*a3*C_A_3*T4a_bar[i]+
			       (Lambda_Q-0.5*a1*C_A_1)*T4b[i]+(Lambda_Q-0.5*a3*C_A_3)*T4b_bar[i]+
			       Lambda_Q*(T7[i]+T7_bar[i]));
	}
      }
    }
    break;
  case UUBG: case DDBG:
     {
        // Finish Group Theory Factors:
        if (process==UUBG) {
           Y_Q = 2.0/3.0;
           Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
           rho13 = Y_Q*sqrt(a1*a3);
        }
        else if (process==DDBG) {
           Y_Q = -1.0/3.0;
           Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
           rho13 = Y_Q*sqrt(a1*a3);
        }
     
        highC.resize(1,5);
     
        for (int i=0; i<5; i++) {
           highC(0,i) = rho13*4.0*pi*(R1[i]+R1_bar[i]);
        
           if (order>=1) {
           
              highC(0,i) += rho13*(-0.5*a3*C_A_3*T1b[i]-0.5*a1*C_A_1*T1b_bar[i]+
                                    (Lambda_Q-0.5*a3*C_A_3-0.5*a1*C_A_1)*(T1c[i]+T1c_bar[i])+
                                    0.5*a3*C_A_3*T3a[i]+0.5*a1*C_A_1*T3a_bar[i]+
                                    (Lambda_Q-0.5*a3*C_A_3)*T3b[i]+(Lambda_Q-0.5*a1*C_A_1)*T3b_bar[i]+
                                    0.5*a1*C_A_1*T4a[i]+0.5*a3*C_A_3*T4a_bar[i]+
                                    (Lambda_Q-0.5*a1*C_A_1)*T4b[i]+(Lambda_Q-0.5*a3*C_A_3)*T4b_bar[i]+
                                    Lambda_Q*(T7[i]+T7_bar[i]));
           }
        }
     }
    break;
  case QQGG: case QtQtGG: case UUGG:
  case tRtRGG: case DDGG:
    {
      // Finish Group Theory Factors:
      if (process==QQGG || process==QtQtGG) {
	Y_Q = 1.0/6.0;
	Lambda_Q = C_F_3*a3 + C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
	Lambda_Q_s = C_F_3*a3;
	Lambda_Q_ew = C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
      }
      else if (process==UUGG || process==tRtRGG) {
	Y_Q = 2.0/3.0;
	Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
	Lambda_Q_s = C_F_3*a3;
	Lambda_Q_ew = Y_Q*Y_Q*C_F_1*a1;
      }
      else if (process==DDGG || process==tRtRGG) {
	Y_Q = -1.0/3.0;
	Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
	Lambda_Q_s = C_F_3*a3;
	Lambda_Q_ew = Y_Q*Y_Q*C_F_1*a1;
      }
      
      highC.resize(3,5);
         
      for (int i=0; i<5; i++) {
	highC(0,i) = G_Plus_SU3*(4.0*pi*a3)*(R1[i]+R1_bar[i]);
	highC(1,i) = G_Plus_SU3_D*(4.0*pi*a3)*(R1[i]+R1_bar[i]);
	highC(2,i) = G_f*4.0*pi*a3*(-0.5*R1[i]+0.5*R1_bar[i]-R2[i]);
        
	if (order>=1) {
	  highC(0,i) += G_Plus_SU3*(a3*(Lambda_Q_ew)*(T1c[i]+T1c_bar[i])+
				    a3*(Lambda_Q_ew)*(T3b[i]+T3b_bar[i])+
				    a3*(Lambda_Q_ew)*(T4b[i]+T4b_bar[i])+
				    a3*Lambda_Q_ew*(T7[i]+T7_bar[i]));
	  highC(1,i) += G_Plus_SU3_D*(a3*(Lambda_Q_ew)*(T1c[i]+T1c_bar[i])+
				      a3*(Lambda_Q_ew)*(T3b[i]+T3b_bar[i])+
				      a3*(Lambda_Q_ew)*(T4b[i]+T4b_bar[i])+
				      a3*Lambda_Q_ew*(T7[i]+T7_bar[i]));
	  highC(2,i) += G_f*(a3*(-0.5*Lambda_Q_ew)*(T1c[i]-T1c_bar[i])-
			     0.5*a3*(Lambda_Q_ew)*(T3b[i]-T3b_bar[i])-
			     0.5*a3*(Lambda_Q_ew)*(T4b[i]-T4b_bar[i])+
			     a3*(Lambda_Q_ew)*T5b[i]-
			     0.5*a3*Lambda_Q_ew*(T7[i]-T7_bar[i]));
	  if (includeAlphaS2) {
	    highC(0,i) += G_Plus_SU3*((-0.5*a3*a3*C_A_3)*(T1b[i]+T1b_bar[i])+
				      a3*(Lambda_Q_s-a3*C_A_3)*(T1c[i]+T1c_bar[i])+
				      0.5*a3*a3*C_A_3*(T3a[i]+T3a_bar[i])+
				      a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T3b[i]+T3b_bar[i])+
				      0.5*a3*a3*C_A_3*(T4a[i]+T4a_bar[i])+
				      a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T4b[i]+T4b_bar[i])+
				      a3*Lambda_Q_s*(T7[i]+T7_bar[i])) + 
	      G_TT*(-a3*a3*(T1a[i]+T1a_bar[i])+a3*a3*(T1b[i]+T1b_bar[i])+
		    a3*a3*(T1c[i]+T1c_bar[i])+2.0*a3*a3*T1d[i]);
	    highC(1,i) += G_Plus_SU3_D*((-0.5*a3*a3*C_A_3)*(T1b[i]+T1b_bar[i])+
					a3*(Lambda_Q_s-a3*C_A_3)*(T1c[i]+T1c_bar[i])+
					0.5*a3*a3*C_A_3*(T3a[i]+T3a_bar[i])+
					a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T3b[i]+T3b_bar[i])+
					0.5*a3*a3*C_A_3*(T4a[i]+T4a_bar[i])+
					a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T4b[i]+T4b_bar[i])+
					a3*Lambda_Q_s*(T7[i]+T7_bar[i])) + 
	      G_TT_3_D*(-a3*a3*(T1a[i]+T1a_bar[i])+a3*a3*(T1b[i]+T1b_bar[i])+
			a3*a3*(T1c[i]+T1c_bar[i])+2.0*a3*a3*T1d[i]);
	    highC(2,i) += G_f*(0.25*a3*a3*C_A_3*(T1a[i]-T1a_bar[i])+
			       a3*(0.25*a3*C_A_3-0.5*Lambda_Q_s)*(T1c[i]-T1c_bar[i])+
			       0.5*a3*a3*C_A_3*T2a[i]+a3*a3*T_CS_3*T2b[i]-
			       0.5*a3*a3*C_A_3*T2c[i]+a3*a3*T_WF_3*T2d[i]-
			       0.25*a3*a3*C_A_3*(T3a[i]-T3a_bar[i])-
			       0.5*a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T3b[i]-T3b_bar[i])-
			       0.25*a3*a3*C_A_3*(T4a[i]-T4a_bar[i])-
			       0.5*a3*(Lambda_Q_s-0.5*a3*C_A_3)*(T4b[i]-T4b_bar[i])+
			       0.5*a3*a3*C_A_3*T5a[i]+a3*(Lambda_Q_s-0.5*a3*C_A_3)*T5b[i]+
			       a3*a3*C_A_3*T6a[i]+a3*a3*C_A_3*T6b[i]+
			       a3*a3*T_CS_3*T6c[i]+a3*a3*T_WF_3*T6d[i]-
			       0.5*a3*Lambda_Q_s*(T7[i]-T7_bar[i]));
	  }
	}
      }
      
      if ( (process==QtQtGG||process==tRtRGG) && order>=1) {
            
	if (process==tRtRGG) {
	  tRorQ = 2.0;
	}
	else {
	  tRorQ = 1.0;
	}
	highC(0,0) += tRorQ*(-1.0*(s*((s+t)*L_t - t*L_u)*y_t*y_t*a3)/(48.*pi*t*u*s));
	highC(0,1) += tRorQ*((s*L_t*y_t*y_t*a3)/(24.*pi*t*s));
	highC(0,2) += tRorQ*(-(s*s*y_t*y_t*a3)/((24.*pi*s*t+24.*pi*t*t)*s));
	highC(0,3) += tRorQ*(-(s*s*y_t*y_t*a3)/((24.*pi*s*t+24.*pi*t*t)*s));
	highC(1,0) += tRorQ*(-1.0*(s*((s+t)*L_t - t*L_u)*y_t*y_t*a3)/(16.*pi*t*u*s));
	highC(1,1) += tRorQ*((s*L_t*y_t*y_t*a3)/(8.*pi*t*s));
	highC(1,2) += tRorQ*(-(s*s*y_t*y_t*a3)/((8.*pi*s*t+8.*pi*t*t)*s));
	highC(1,3) += tRorQ*(-(s*s*y_t*y_t*a3)/((8.*pi*s*t+8.*pi*t*t)*s));
	highC(2,0) += tRorQ*((s*((s+t)*L_t + t*L_u)*y_t*y_t*a3)/(16.*pi*t*u*s));
	highC(2,1) += tRorQ*(((2.*t-2.*t*L_s-s*L_t)*y_t*y_t*a3)/(8.*pi*t*s));
	highC(2,2) += tRorQ*(-1.0*(s*(s+2.*t)*y_t*y_t*a3)/(8.*pi*t*u*s));
	highC(2,3) += tRorQ*(-1.0*(s*(s+2.*t)*y_t*y_t*a3)/(8.*pi*t*u*s));
      }
    }
    break;
  default:
    assert(false);
  }  
  return highC;
}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
HighEnergyMatching::Spin0HighMatching(Energy highScale, 
				      Energy2 s, Energy2 t, Energy2 u,
				      EWProcess::Process process, 
				      bool oneLoop, bool ) {
  using Constants::pi;
  unsigned int order = !oneLoop? 0 : 1;
  // (If crossed graphs, swap s and t here)
  Complex L_s = MinusLog(-s/(highScale*highScale));
  Complex L_t = MinusLog(-t/(highScale*highScale));
  Complex L_u = MinusLog(-u/(highScale*highScale));
  Complex L_s2 = L_s*L_s;
  Complex L_t2 = L_t*L_t;
  Complex L_u2 = L_u*L_u;

  // Tree-Level:
  complex<InvEnergy2> S1 = 2.0/s;
   
  // Topology T1:
  complex<InvEnergy2> T1b = (-L_s2/(2.0*u)*(7.0*t/s+3.0)+2.0/u*L_t2+L_s*L_t*4.0/u*(t-u)/s+
			     L_s*2.0/s-4.0/s-pi*pi/(4.0*u)*(11.0+19.0*t/s));
  complex<InvEnergy2> T1b_bar = -1.0*(-L_s2/(2.0*t)*(7.0*u/s+3.0)+2.0/t*L_u2+L_s*L_u*4.0/t*(u-t)/s+
				      L_s*2.0/s-4.0/s-pi*pi/(4.0*t)*(11.0+19.0*u/s));
   
  // Topologies T2:
  complex<InvEnergy2> T2a = 1.0/s*(-2.0*L_s2+8.0*L_s-16.0+pi*pi/3.0);
  complex<InvEnergy2> T2b = 1.0/s*(0.5*L_s2+2.0*L_s-4.0-pi*pi/12.0);
  
  // Topologies T5:
  complex<InvEnergy2> T5a = 1.0/s*(-2.0*L_s2+6.0*L_s+pi*pi/3.0-16.0);
  complex<InvEnergy2> T5b = 1.0/s*(2.0*L_s-4.0);
  
  // Topologies T6:
  complex<InvEnergy2> T6a = 1.0/s*(-19.0/6.0*L_s+58.0/9.0);
  complex<InvEnergy2> T6b = 1.0/s*(-1.0/6.0*L_s+4.0/9.0);
  complex<InvEnergy2> T6c = 1.0/s*(2.0/3.0*L_s-16.0/9.0);
  complex<InvEnergy2> T6d = 1.0/s*(4.0/3.0*L_s-20.0/9.0);
   
  // Group Theory Factors / SM parameters needed for matrix elements:
  double a1 = ElectroWeakReweighter::coupling()->a1(highScale);
  double a2 = ElectroWeakReweighter::coupling()->a2(highScale);
  double a3 = ElectroWeakReweighter::coupling()->a3(highScale);
  double y_t = ElectroWeakReweighter::coupling()->y_t(highScale);
  double Y_phi = 1.0/2.0;
  double C_F_3 = 4.0/3.0;
  double C_F_2 = 3.0/4.0;
  double C_F_1 = 1.0;
  double n_g = 3.0;
  double n_S = 1.0;
  // Factors TBD after fermion helicity is specified:
  double Y_Q(0.), Lambda_Q(0.), Lambda_phi(0.);
  boost::numeric::ublas::matrix<complex<InvEnergy2> > highC(1,1);
  switch (process) {
    
  case QQPhiPhi: case LLPhiPhi:
    // Finish Group Theory Factors:
    if (process==QQPhiPhi) {
      Y_Q = 1.0/6.0;
      Lambda_Q = C_F_3*a3 + C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
      Lambda_phi = C_F_2*a2+Y_phi*Y_phi*a1;
    }
    else if (process==LLPhiPhi) {
      Y_Q = -1.0/2.0;
      Lambda_Q = C_F_2*a2 + Y_Q*Y_Q*C_F_1*a1;
      Lambda_phi = C_F_2*a2+Y_phi*Y_phi*a1;
    }
    highC.resize(2,1);
    highC(0,0) = S1*(4.0*pi*a2);
    highC(1,0) = S1*(4.0*pi*a1*Y_Q*Y_phi);
    if (order>=1) {
      highC(0,0) += T1b*(0.5*a2*a2+2.0*a1*a2*Y_Q*Y_phi) + 
	T1b_bar*(-0.5*a2*a2+2.0*a1*a2*Y_Q*Y_phi) + 
	T2a*(-a2*a2+Lambda_phi*a2) + T2b*a2*a2 + 
	T5a*(-a2*a2+Lambda_Q*a2) + T5b*a2*a2 + 
	T6a*2.0*a2*a2 + T6b*2.0*a2*a2 + 
	T6c*0.5*a2*a2*n_S + T6d*2.0*a2*a2*n_g;
      highC(1,0) += T1b*(3.0/16.0*a2*a2+a1*a1*Y_Q*Y_Q*Y_phi*Y_phi) + 
	T1b_bar*(3.0/16.0*a2*a2+a1*a1*Y_Q*Y_Q*Y_phi*Y_phi) + 
	T2a*(Lambda_phi*a1*Y_Q*Y_phi) + T5a*(Lambda_Q*a1*Y_Q*Y_phi) + 
	T6c*(2.0*a1*a1*n_S*Y_Q*Y_phi*Y_phi*Y_phi) + 
	T6d*(10.0/3.0*a1*a1*n_g*Y_Q*Y_phi);
      // Top Quark contributions:
      highC(0,0) += -3.0*y_t*y_t*a2/(4.0*pi)/s*(2.0*L_s-4.0);
      highC(1,0) += -3.0*y_t*y_t*a1/(4.0*pi)*(Y_Q*Y_phi)/s*(2.0*L_s-4.0);
    }
    break;
  case UUPhiPhi: case DDPhiPhi: case EEPhiPhi:
    // Finish Group Theory Factors:
    if (process==UUPhiPhi) {
      Y_Q = 2.0/3.0;
      Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
      Lambda_phi = C_F_2*a2 + Y_phi*Y_phi*a1;
    }
    else if (process==DDPhiPhi) {
      Y_Q = -1.0/3.0;
      Lambda_Q = C_F_3*a3 + Y_Q*Y_Q*C_F_1*a1;
      Lambda_phi = C_F_2*a2 + Y_phi*Y_phi*a1;
    }
    else if (process==EEPhiPhi) {
      Y_Q = -1.0;
      Lambda_Q = Y_Q*Y_Q*C_F_1*a1;
      Lambda_phi = C_F_2*a2 + Y_phi*Y_phi*a1;
    }
    
    highC.resize(1,1);
    highC(0,0) = ZERO;
    
    highC(0,0) = S1*(4.0*pi*a1*Y_Q*Y_phi);
    
    if (order>=1) {
      highC(0,0) += T1b*(a1*a1*Y_Q*Y_Q*Y_phi*Y_phi) + 
	T1b_bar*(a1*a1*Y_Q*Y_Q*Y_phi*Y_phi) + 
	T2a*(Lambda_phi*a1*Y_Q*Y_phi) + T5a*(Lambda_Q*a1*Y_Q*Y_phi) + 
	T6c*(2.0*a1*a1*n_S*Y_Q*Y_phi*Y_phi*Y_phi) + 
	T6d*(10.0/3.0*a1*a1*n_g*Y_Q*Y_phi);
      // Top Quark Contribution:
      highC(0,0) += -3.0*y_t*y_t*a1/(4.0*pi)*(Y_Q*Y_phi)/s*(2.0*L_s-4.0);
    }
    break;
  default:
    assert(false);
  }
  return highC;
}
