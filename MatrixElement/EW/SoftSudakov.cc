// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SoftSudakov class.
//

#include "SoftSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "GroupInvariants.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace GroupInvariants;

SoftSudakov::SoftSudakov() : K_ORDER_(3) {}

SoftSudakov::~SoftSudakov() {}

IBPtr SoftSudakov::clone() const {
  return new_ptr(*this);
}

IBPtr SoftSudakov::fullclone() const {
  return new_ptr(*this);
}

void SoftSudakov::persistentOutput(PersistentOStream & os) const {
  os << K_ORDER_;
}

void SoftSudakov::persistentInput(PersistentIStream & is, int) {
  is >> K_ORDER_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SoftSudakov,Interfaced>
describeHerwigSoftSudakov("Herwig::SoftSudakov", "HwMEEW.so");

void SoftSudakov::Init() {

  static ClassDocumentation<SoftSudakov> documentation
    ("The SoftSudakov class implements the soft EW Sudakov");

}

InvEnergy SoftSudakov::operator ()(Energy mu) const {
  // Include K-factor Contributions (Cusps):
  GaugeContributions cusp = cuspContributions(mu,K_ORDER_,high_);
  Complex gamma = cusp.SU3*G3_(row_,col_) + cusp.SU2*G2_(row_,col_) + cusp.U1*G1_(row_,col_);
  if (real_) {
    return gamma.real()/mu;
  }
  else {
    return gamma.imag()/mu;
  }
}

boost::numeric::ublas::matrix<Complex> 
SoftSudakov::evaluateSoft(boost::numeric::ublas::matrix<Complex> & G3,
			  boost::numeric::ublas::matrix<Complex> & G2,
			  boost::numeric::ublas::matrix<Complex> & G1,
			  Energy mu_h, Energy mu_l, bool high) {
  assert( G3.size1() == G2.size1() && G3.size1() == G1.size1() &&
	  G3.size2() == G2.size2() && G3.size2() == G1.size2() && 
	  G3.size1() == G3.size2());
  G3_ = G3;
  G2_ = G2;
  G1_ = G1;
  high_ = high;
  unsigned int NN = G3_.size1();
  // gamma is the matrix to be numerically integrated to run the coefficients.
  boost::numeric::ublas::matrix<Complex> gamma(NN,NN);
  for(row_=0;row_<NN;++row_) {
    for(col_=0;col_<NN;++col_) {
      real_ = true;
      gamma(row_,col_).real(integrator_.value(*this,mu_h,mu_l));
      real_ = false;
      gamma(row_,col_).imag(integrator_.value(*this,mu_h,mu_l));
    }
  }
  // Resummed:
  //return gamma.exp();
  // Fixed Order:
  return boost::numeric::ublas::identity_matrix<Complex>(NN,NN) + gamma;
}

boost::numeric::ublas::matrix<Complex>
SoftSudakov::lowEnergyRunning(Energy EWScale, Energy lowScale, 
			      Energy2 s, Energy2 t, Energy2 u, 
			      Herwig::EWProcess::Process process) {
  using namespace EWProcess;
  using Constants::pi;
  static const Complex I(0,1.0);
  Complex T = getT(s,t), U = getU(s,u);
  boost::numeric::ublas::matrix<Complex> G1, G2, G3;
  unsigned int numBrokenGauge;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    {
      numBrokenGauge = 12;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      boost::numeric::ublas::matrix<Complex> gam3 = Gamma3(U,T);
      for (unsigned int i=0; i<numBrokenGauge/2; i++) {
	G3(i,i)     += gam3(0,0);
	G3(i,i+6)   += gam3(0,1);
	G3(i+6,i)   += gam3(1,0);
	G3(i+6,i+6) += gam3(1,1);
      }
      G1(0,0) = G1(6,6)   = Gamma1(2.0/3.0,2.0/3.0,2.0/3.0,2.0/3.0,T,U);
      G1(1,1) = G1(7,7)   = Gamma1(-1.0/3.0,-1.0/3.0,2.0/3.0,2.0/3.0,T,U);
      G1(2,2) = G1(8,8)   = Gamma1(2.0/3.0,2.0/3.0,-1.0/3.0,-1.0/3.0,T,U);
      G1(3,3) = G1(9,9)   = Gamma1(-1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,T,U);
      G1(4,4) = G1(10,10) = Gamma1(-1.0/3.0,2.0/3.0,2.0/3.0,-1.0/3.0,T,U);
      G1(5,5) = G1(11,11) = Gamma1(2.0/3.0,-1.0/3.0,-1.0/3.0,2.0/3.0,T,U);
    }
    break;
  case QQUU:
  case QtQtUU:
  case QQtRtR:
    {
      numBrokenGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      boost::numeric::ublas::matrix<Complex> gam3 = Gamma3(U,T);
      for (unsigned int i=0; i<numBrokenGauge/2; i++) {
	G3(i,i)     += gam3(0,0);
	G3(i,i+2)   += gam3(0,1);
	G3(i+2,i)   += gam3(1,0);
	G3(i+2,i+2) += gam3(1,1);
      }
      G1(0,0) = G1(2,2) = Gamma1(2.0/3.0,2.0/3.0,T,U);
      G1(1,1) = G1(3,3) = Gamma1(2.0/3.0,-1.0/3.0,T,U);
    }
    break;
  case QQDD:
  case QtQtDD:
    {
      numBrokenGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      boost::numeric::ublas::matrix<Complex> gam3 = Gamma3(U,T);
      for (unsigned int i=0; i<numBrokenGauge/2; i++) {
	G3(i,i) += gam3(0,0);
	G3(i,i+2) += gam3(0,1);
	G3(i+2,i) += gam3(1,0);
	G3(i+2,i+2) += gam3(1,1);
      }
      G1(0,0) = G1(2,2) = Gamma1(-1.0/3.0,2.0/3.0,T,U);
      G1(1,1) = G1(3,3) = Gamma1(-1.0/3.0,-1.0/3.0,T,U);
    }
    break;
  case QQLL:
    {
      numBrokenGauge = 6;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
            G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(2.0/3.0,2.0/3.0,0.0,0.0,T,U);
      G1(1,1) = Gamma1(-1.0/3.0,-1.0/3.0,0.0,0.0,T,U);
      G1(2,2) = Gamma1(2.0/3.0,2.0/3.0,-1.0,-1.0,T,U);
      G1(3,3) = Gamma1(-1.0/3.0,-1.0/3.0,-1.0,-1.0,T,U);
      G1(4,4) = Gamma1(-1.0/3.0,2.0/3.0,0.0,-1.0,T,U);
      G1(5,5) = Gamma1(2.0/3.0,-1.0/3.0,-1.0,0.0,T,U);
    }
    break;
  case QQEE:
    {
      numBrokenGauge = 2;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G1 *= 0.0; G2 *= 0.0; G3 *= 0.0;
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<2; i++) {
	G3(i,i) += gam3s;
      }
      G1(0,0) = Gamma1(2.0/3.0,-1.0,T,U);
      G1(1,1) = Gamma1(-1.0/3.0,-1.0,T,U);
    }
    break;
  case UUUU:
  case UUUUiden:
  case tRtRUU:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1 *= 0.0; G2 *= 0.0; G3 *= 0.0;
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,2.0/3.0,T,U);
    break;
  case UUDD:
  case tRtRDD:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,2.0/3.0,T,U);
    break;
  case UULL:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = G3(1,1) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0,0.0,T,U);
    G1(1,1) = Gamma1(2.0/3.0,-1.0,T,U);
    break;
  case UUEE:
    numBrokenGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0,-1.0,T,U);
    break;
  case DDDD:
  case DDDDiden:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,-1.0/3.0,T,U);
    break;
  case DDLL:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = G3(1,1) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0,0.0,T,U);
    G1(1,1) = Gamma1(-1.0/3.0,-1.0,T,U);
    break;
  case DDEE:
    numBrokenGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0,-1.0,T,U);
    break;
  case LLLL:
  case LLLLiden:
    numBrokenGauge = 6;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(0.0,0.0,0.0,0.0,T,U);
    G1(1,1) = Gamma1(-1.0,-1.0,0.0,0.0,T,U);
    G1(2,2) = Gamma1(0.0,0.0,-1.0,-1.0,T,U);
    G1(3,3) = Gamma1(-1.0,-1.0,-1.0,-1.0,T,U);
    G1(4,4) = Gamma1(-1.0,0.0,0.0,-1.0,T,U);
    G1(5,5) = Gamma1(0.0,-1.0,-1.0,0.0,T,U);
    break;
    
  case LLEE:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(0.0,-1.0,T,U);
    G1(1,1) = Gamma1(-1.0,-1.0,T,U);
    break;
  case EEEE:
  case EEEEiden:
    numBrokenGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(-1.0,-1.0,T,U);
    break;
  case QQWW:
    {
      numBrokenGauge = 20;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(2./3.,2./3.,-1.,-1.,T,U);
      G1(1,1) = Gamma1(2./3.,2./3.,1.,1.,T,U);
      G1(2,2) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(3,3) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(4,4) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(5,5) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(6,6) = Gamma1(-1./3.,-1./3.,-1.,-1.,T,U);
      G1(7,7) = Gamma1(-1./3.,-1./3.,1.,1.,T,U);
      G1(8,8) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(9,9) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(10,10) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(11,11) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(12,12) = Gamma1(-1./3.,2./3.,0.,-1.,T,U);
      G1(13,13) = Gamma1(-1./3.,2./3.,0.,-1.,T,U);
      G1(14,14) = Gamma1(-1./3.,2./3.,1.,0.,T,U);
      G1(15,15) = Gamma1(-1./3.,2./3.,1.,0.,T,U);
      G1(16,16) = Gamma1(2./3.,-1./3.,-1.,0.,T,U);
      G1(17,17) = Gamma1(2./3.,-1./3.,-1.,0.,T,U);
      G1(18,18) = Gamma1(2./3.,-1./3.,0.,1.,T,U);
      G1(19,19) = Gamma1(2./3.,-1./3.,0.,1.,T,U);
    }
    break;
  case QQPhiPhi:
    {
      numBrokenGauge = 14;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(2./3.,2./3.,1.,1.,T,U);
      G1(1,1) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(2,2) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(3,3) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(4,4) = Gamma1(2./3.,2./3.,0.,0.,T,U);
      G1(5,5) = Gamma1(-1./3.,-1./3.,1.,1.,T,U);
      G1(6,6) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(7,7) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(8,8) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(9,9) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
      G1(10,10) = Gamma1(-1./3.,2./3.,1.,0.,T,U);
      G1(11,11) = Gamma1(-1./3.,2./3.,1.,0.,T,U);
      G1(12,12) = Gamma1(2./3.,-1./3.,0.,1.,T,U);
      G1(13,13) = Gamma1(2./3.,-1./3.,0.,1.,T,U);
    }
    break;
  case QQWG:
    numBrokenGauge = 6;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    for (unsigned int i=0; i<numBrokenGauge; i++) {
      G3(i,i) = -17.0/6.0*I*pi + 3.0/2.0*(U+T);
    }
    G1(0,0) = Gamma1(-1./3.,2./3.,1.,0.,T,U);
    G1(1,1) = Gamma1(2./3.,-1./3.,-1.,0.,T,U);
    G1(2,2) = Gamma1(2./3.,2./3.,0.,0.,T,U);
    G1(3,3) = Gamma1(2./3.,2./3.,0.,0.,T,U);
    G1(4,4) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
    G1(5,5) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
    break;
  case QQBG:
    numBrokenGauge = 4;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    for (unsigned int i=0; i<numBrokenGauge; i++) {
      G3(i,i) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    }
    G1(0,0) = Gamma1(2./3.,2./3.,0.,0.,T,U);
    G1(1,1) = Gamma1(2./3.,2./3.,0.,0.,T,U);
    G1(2,2) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
    G1(3,3) = Gamma1(-1./3.,-1./3.,0.,0.,T,U);
    break;
  case QQGG:
  case QtQtGG:
    {
      numBrokenGauge = 6;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      boost::numeric::ublas::matrix<Complex> gam3g = Gamma3g(U,T);
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  G3(ix  ,iy  ) = gam3g(ix,iy);
	  G3(ix+3,iy+3) = gam3g(ix,iy);
	}
      }
      G1(0,0) = G1(1,1) = G1(2,2) = Gamma1(2./3.,0.,T,U);
      G1(3,3) = G1(4,4) = G1(5,5) = Gamma1(-1./3.,0.,T,U);
    }
    break;
  case LLWW:
    numBrokenGauge = 20;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(0.,0.,-1.,-1.,T,U);
    G1(1,1) = Gamma1(0.,0.,1.,1.,T,U);
    G1(2,2) = Gamma1(0.,0.,0.,0.,T,U);
    G1(3,3) = Gamma1(0.,0.,0.,0.,T,U);
    G1(4,4) = Gamma1(0.,0.,0.,0.,T,U);
    G1(5,5) = Gamma1(0.,0.,0.,0.,T,U);
    G1(6,6) = Gamma1(-1.,-1.,-1.,-1.,T,U);
    G1(7,7) = Gamma1(-1.,-1.,1.,1.,T,U);
    G1(8,8) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(9,9) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(10,10) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(11,11) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(12,12) = Gamma1(-1.,0.,0.,-1.,T,U);
    G1(13,13) = Gamma1(-1.,0.,0.,-1.,T,U);
    G1(14,14) = Gamma1(-1.,0.,1.,0.,T,U);
    G1(15,15) = Gamma1(-1.,0.,1.,0.,T,U);
    G1(16,16) = Gamma1(0.,-1.,-1.,0.,T,U);
    G1(17,17) = Gamma1(0.,-1.,-1.,0.,T,U);
    G1(18,18) = Gamma1(0.,-1.,0.,1.,T,U);
    G1(19,19) = Gamma1(0.,-1.,0.,1.,T,U);
    break;
  case LLPhiPhi:
    numBrokenGauge = 14;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1 *= 0.0; G2 *= 0.0; G3 *= 0.0;
    G1(0,0) = Gamma1(0.,0.,1.,1.,T,U);
    G1(1,1) = Gamma1(0.,0.,0.,0.,T,U);
    G1(2,2) = Gamma1(0.,0.,0.,0.,T,U);
    G1(3,3) = Gamma1(0.,0.,0.,0.,T,U);
    G1(4,4) = Gamma1(0.,0.,0.,0.,T,U);
    G1(5,5) = Gamma1(-1.,-1.,1.,1.,T,U);
    G1(6,6) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(7,7) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(8,8) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(9,9) = Gamma1(-1.,-1.,0.,0.,T,U);
    G1(10,10) = Gamma1(-1.,0.,1.,0.,T,U);
    G1(11,11) = Gamma1(-1.,0.,1.,0.,T,U);
    G1(12,12) = Gamma1(0.,-1.,0.,1.,T,U);
    G1(13,13) = Gamma1(0.,-1.,0.,1.,T,U);
    break;
  case UUBB:
    {
      numBrokenGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
	G1(i,i) = Gamma1(2./3.,0.,T,U);
      }
    }
    break;
  case UUPhiPhi:
    {
      numBrokenGauge = 5;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(2.0/3.0,2.0/3.0,1.,1.,T,U);
      G1(1,1) = G1(2,2) = G1(3,3) = G1(4,4) = Gamma1(2./3.,0.,T,U);
    }
    break;
  case UUBG:
    {
      numBrokenGauge = 2;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam1 = Gamma1(2./3.,0.,T,U);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
	G1(i,i) = gam1;
      }
    }
    break;
  case UUGG:
  case tRtRGG:
    numBrokenGauge = 3;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3g(U,T);
    G1(0,0) = G1(1,1) = G1(2,2) = Gamma1(2./3.,0.,T,U);
    break;
  case DDBB:
    {
      numBrokenGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      Complex gam1  = Gamma1(-1./3.,0.,T,U);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
	G1(i,i) = gam1;
      }
    }
    break;
  case DDPhiPhi:
    {
      numBrokenGauge = 5;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(-1.0/3.0,-1.0/3.0,1.,1.,T,U);
      G1(1,1) = G1(2,2) = G1(3,3) = G1(4,4) = Gamma1(-1./3.,0.,T,U);
    }
    break;
  case DDBG:
    numBrokenGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    for (unsigned int i=0; i<numBrokenGauge; i++) {
      G3(i,i) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    }
    G1(0,0) = G1(1,1) = Gamma1(-1./3.,0.,T,U);
    break;
  case DDGG:
    numBrokenGauge = 3;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3g(U,T);
    G1(0,0) = G1(1,1) = G1(2,2) = Gamma1(-1./3.,0.,T,U);
    break;
  case EEBB:
    {
      numBrokenGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam1 = Gamma1(-1.,0.,T,U);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G1(i,i) = gam1;
      };
    }
    break;
  case EEPhiPhi:
    numBrokenGauge = 5;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(-1.,1.,T,U);
    G1(1,1) = G1(2,2) = G1(3,3) = G1(4,4) = Gamma1(-1.,0.,T,U);
    break;         
  default:
    assert(false);
  }
   
  if (EWScale==lowScale) {
    return boost::numeric::ublas::identity_matrix<Complex>(G1.size1());
  }
  else {
    return evaluateSoft(G3,G2,G1,EWScale,lowScale,false);
  }
}

boost::numeric::ublas::matrix<Complex>
SoftSudakov::highEnergyRunning(Energy highScale, Energy EWScale,
			       Energy2 s, Energy2 t, Energy2 u,
			       Herwig::EWProcess::Process process) {
  using namespace EWProcess;
  using Constants::pi;
  static const Complex I(0,1.0);
  Complex T = getT(s,t), U = getU(s,u);
  boost::numeric::ublas::matrix<Complex> G1,G2,G3;
  unsigned int numGauge;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    {
      numGauge = 4;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      boost::numeric::ublas::matrix<Complex> gam3 = Gamma3(U,T);
      G3(0,0) += gam3(0,0);
      G3(0,2) += gam3(0,1);
      G3(2,0) += gam3(1,0);
      G3(2,2) += gam3(1,1);
      G3(1,1) += gam3(0,0);
      G3(1,3) += gam3(0,1);
      G3(3,1) += gam3(1,0);
      G3(3,3) += gam3(1,1);
      boost::numeric::ublas::matrix<Complex> gam2 = Gamma2(U,T);
      G2(0,0) += gam2(0,0);
      G2(0,1) += gam2(0,1);
      G2(1,0) += gam2(1,0);
      G2(1,1) += gam2(1,1);
      G2(2,2) += gam2(0,0);
      G2(2,3) += gam2(0,1);
      G2(3,2) += gam2(1,0);
      G2(3,3) += gam2(1,1);
      G1(0,0) = G1(1,1) = G1(2,2) = G1(3,3) = Gamma1(1.0/6.0,1.0/6.0,T,U);
    }
    break;
  case QQUU:
  case QtQtUU:
  case QQtRtR:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G2 = Gamma2Singlet();
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,1.0/6.0,T,U);
    break;
  case QQDD:
  case QtQtDD:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G2 = Gamma2Singlet();
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,1.0/6.0,T,U);
    break;
  case QQLL:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3Singlet();
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/2.0,1.0/6.0,T,U);
    break;
  case QQEE:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge); 
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,1.0/6.0,T,U);
    break;
  case UUUU:
  case UUUUiden:
  case tRtRUU:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,2.0/3.0,T,U);
    break;
  case UUDD:
  case tRtRDD:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,2.0/3.0,T,U);
    break;
  case UULL:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/2.0,2.0/3.0,T,U);
    break;
  case UUEE:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,2.0/3.0,T,U);
    break;
  case DDDD:
  case DDDDiden:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,-1.0/3.0,T,U);
    break;
  case DDLL:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/2.0,-1.0/3.0,T,U);
    break;
  case DDEE:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,-1.0/3.0,T,U);
    break;
  case LLLL:
  case LLLLiden:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/2.0,-1.0/2.0,T,U);
    break;
  case LLEE:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,-1.0/2.0,T,U);
    break;
  case EEEE:
  case EEEEiden:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G1(0,0) = Gamma1(-1.0,-1.0,T,U);
    break;
  case QQWW:
    {
      numGauge = 5;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<5; i++) {
	G3(i,i) = gam3s;
	G1(i,i) = Gamma1(1.0/6.0);
      }
      G2 = Gamma2w(U,T);
    }
    break;
  case QQPhiPhi:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3Singlet();
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(1.0/2.0,1.0/6.0,T,U);
    break;
  case QQWG:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -17.0/6.0*I*pi + 3.0/2.0*(U+T);
    G2(0,0) = -7.0/4.0*I*pi + (U+T);
    G1(0,0) = Gamma1(1.0/6.0);
    break;
  case QQBG:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G2(0,0) = -3.0/4.0*I*pi;
    G1(0,0) = Gamma1(1.0/6.0);
    break;
  case QQGG:
  case QtQtGG:
    {
      numGauge = 3;
      G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
      G3 = Gamma3g(T,U);
      Complex gam2s = Gamma2Singlet()(0,0);
      for (unsigned int i=0; i<3; i++) {
	G2(i,i) = gam2s;
	G1(i,i) = Gamma1(1.0/6.0);
      }
    }
    break;
  case LLWW:
    numGauge = 5;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    for (unsigned int i=0; i<5; i++) {
      G1(i,i) = Gamma1(-1.0/2.0);
    }
    G2 = Gamma2w(U,T);
    break;
  case LLPhiPhi:
    numGauge = 2;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(1.0/2.0,-1.0/2.0,T,U);
    break;
  case UUBB:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0);
    break;
  case UUPhiPhi:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,2.0/3.0,T,U);
    break;
  case UUBG:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G1(0,0) = Gamma1(2.0/3.0);
    break;
  case UUGG:
  case tRtRGG:
    numGauge = 3;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3g(U,T);
    for (unsigned int i=0; i<3; i++) {
      G1(i,i) = Gamma1(2.0/3.0);
    }
    break;
  case DDBB:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0);
    break;
  case DDPhiPhi:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,-1.0/3.0,T,U);
    break;
  case DDBG:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G1(0,0) = Gamma1(-1.0/3.0);
    break;
  case DDGG:
    numGauge = 3;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge); 
    G3 = Gamma3g(U,T);
    for (unsigned int i=0; i<3; i++) {
      G1(i,i) = Gamma1(-1.0/3.0);
    }
    break;
  case EEBB:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G1(0,0) = Gamma1(-1.0);
    break;
  case EEPhiPhi:
    numGauge = 1;
    G1 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G3 = boost::numeric::ublas::zero_matrix<Complex>(numGauge,numGauge);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,-1.0,T,U);
    break;
  default:
    assert(false);
  }
  return evaluateSoft(G3,G2,G1,highScale,EWScale,true);
}
