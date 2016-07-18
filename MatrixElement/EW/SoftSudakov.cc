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
#include "expm-1.h"

using namespace Herwig;
using namespace GroupInvariants;

SoftSudakov::SoftSudakov() : K_ORDER_(3), integrator_(0.,1e-5,1000) {}

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
      if(G3_(row_,col_) == 0. && G2_(row_,col_) == 0. && G1_(row_,col_) == 0.) {
	gamma(row_,col_) = 0.;
      }
      else {
	real_ = true;
	gamma(row_,col_).real(integrator_.value(*this,mu_h,mu_l));
	real_ = false;
	gamma(row_,col_).imag(integrator_.value(*this,mu_h,mu_l));
      }
    }
  }
  // Resummed:
  return boost::numeric::ublas::expm_pad(gamma,7);
}

boost::numeric::ublas::matrix<Complex>
SoftSudakov::lowEnergyRunning(Energy EWScale, Energy lowScale, 
			      Energy2 s, Energy2 t, Energy2 u, 
			      Herwig::EWProcess::Process process,
			      unsigned int iswap) {
  using namespace EWProcess;
  using namespace boost::numeric::ublas;
  using Constants::pi;
  static const Complex I(0,1.0);
  Complex T = getT(s,t), U = getU(s,u);
  matrix<Complex> G1, G2, G3;
  unsigned int numBrokenGauge;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    {
      assert(iswap==0);
      numBrokenGauge = 12;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      matrix<Complex> gam3 = Gamma3(U,T);
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
      assert(iswap==0);
      numBrokenGauge = 4;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      matrix<Complex> gam3 = Gamma3(U,T);
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
      assert(iswap==0);
      numBrokenGauge = 4;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      matrix<Complex> gam3 = Gamma3(U,T);
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
      assert(iswap==0);
      numBrokenGauge = 6;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      assert(iswap==0);
      numBrokenGauge = 2;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,2.0/3.0,T,U);
    break;
  case UUDD:
  case tRtRDD:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,2.0/3.0,T,U);
    break;
  case UULL:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = G3(1,1) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0,0.0,T,U);
    G1(1,1) = Gamma1(2.0/3.0,-1.0,T,U);
    break;
  case UUEE:
    assert(iswap==0);
    numBrokenGauge = 1;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0,-1.0,T,U);
    break;
  case DDDD:
  case DDDDiden:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,-1.0/3.0,T,U);
    break;
  case DDLL:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = G3(1,1) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0,0.0,T,U);
    G1(1,1) = Gamma1(-1.0/3.0,-1.0,T,U);
    break;
  case DDEE:
    assert(iswap==0);
    numBrokenGauge = 1;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0,-1.0,T,U);
    break;
  case LLLL:
  case LLLLiden:
    assert(iswap==0);
    numBrokenGauge = 6;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(0.0,0.0,0.0,0.0,T,U);
    G1(1,1) = Gamma1(-1.0,-1.0,0.0,0.0,T,U);
    G1(2,2) = Gamma1(0.0,0.0,-1.0,-1.0,T,U);
    G1(3,3) = Gamma1(-1.0,-1.0,-1.0,-1.0,T,U);
    G1(4,4) = Gamma1(-1.0,0.0,0.0,-1.0,T,U);
    G1(5,5) = Gamma1(0.0,-1.0,-1.0,0.0,T,U);
    break;
    
  case LLEE:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(0.0,-1.0,T,U);
    G1(1,1) = Gamma1(-1.0,-1.0,T,U);
    break;
  case EEEE:
  case EEEEiden:
    assert(iswap==0);
    numBrokenGauge = 1;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(-1.0,-1.0,T,U);
    break;
  case QQWW:
    {
      assert(iswap==0);
      numBrokenGauge = 20;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      assert(iswap==0);
      numBrokenGauge = 14;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
    assert(iswap==0);
    numBrokenGauge = 6;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
    assert(iswap==0);
    numBrokenGauge = 4;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      matrix<Complex> gam3g;
      Complex gam1a(0.),gam1b(0.);
      if(iswap==0) {
	gam3g = Gamma3g(U,T);
	gam1a = Gamma1( 2./3.,0.,T,U);
	gam1b = Gamma1(-1./3.,0.,T,U);
      }
      else if(iswap==1) {
	gam3g = Gamma3gST(U,T);
	gam1a = Gamma1ST( 2./3.,0.,T,U);
	gam1b = Gamma1ST(-1./3.,0.,T,U);
      }
      else
	assert(false);
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  G3(ix  ,iy  ) = gam3g(ix,iy);
	  G3(ix+3,iy+3) = gam3g(ix,iy);
	}
      }
      G1(0,0) = G1(1,1) = G1(2,2) = gam1a;
      G1(3,3) = G1(4,4) = G1(5,5) = gam1b;
    }
    break;
  case LLWW:
    assert(iswap==0);
    numBrokenGauge = 20;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
    assert(iswap==0);
    numBrokenGauge = 14;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      assert(iswap==0);
      numBrokenGauge = 4;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
	G1(i,i) = Gamma1(2./3.,0.,T,U);
      }
    }
    break;
  case UUPhiPhi:
    {
      assert(iswap==0);
      numBrokenGauge = 5;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      assert(iswap==0);
      numBrokenGauge = 2;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam1 = Gamma1(2./3.,0.,T,U);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
	G1(i,i) = gam1;
      }
    }
    break;
  case DDGG:
  case UUGG:
  case tRtRGG:
    {
      numBrokenGauge = 3;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam1(0.);
      double Y = process==DDGG ? -1./3. : 2./3.;
      if(iswap==0) {
	G3   = Gamma3g(U,T);
	gam1 = Gamma1(Y,0.,T,U);
      }
      else if(iswap==1) {
	G3   = Gamma3gST(U,T);
	gam1 = Gamma1ST(Y,0.,T,U);
      }
      else
	assert(false);    
      G1(0,0) = G1(1,1) = G1(2,2) = gam1;
    }
    break;
  case DDBB:
    {
      assert(iswap==0);
      numBrokenGauge = 4;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
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
      assert(iswap==0);
      numBrokenGauge = 5;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G3(i,i) = gam3s;
      }
      G1(0,0) = Gamma1(-1.0/3.0,-1.0/3.0,1.,1.,T,U);
      G1(1,1) = G1(2,2) = G1(3,3) = G1(4,4) = Gamma1(-1./3.,0.,T,U);
    }
    break;
  case DDBG:
    assert(iswap==0);
    numBrokenGauge = 2;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    for (unsigned int i=0; i<numBrokenGauge; i++) {
      G3(i,i) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    }
    G1(0,0) = G1(1,1) = Gamma1(-1./3.,0.,T,U);
    break;
  case EEBB:
    {
      assert(iswap==0);
      numBrokenGauge = 4;
      G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
      Complex gam1 = Gamma1(-1.,0.,T,U);
      for (unsigned int i=0; i<numBrokenGauge; i++) {
	G1(i,i) = gam1;
      };
    }
    break;
  case EEPhiPhi:
    assert(iswap==0);
    numBrokenGauge = 5;
    G1 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G2 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G3 = zero_matrix<Complex>(numBrokenGauge,numBrokenGauge);
    G1(0,0) = Gamma1(-1.,1.,T,U);
    G1(1,1) = G1(2,2) = G1(3,3) = G1(4,4) = Gamma1(-1.,0.,T,U);
    break;         
  default:
    assert(false);
  }
  // return the answer
  if (EWScale==lowScale) {
    return identity_matrix<Complex>(G1.size1());
  }
  else {
    return evaluateSoft(G3,G2,G1,EWScale,lowScale,false);
  }
}

boost::numeric::ublas::matrix<Complex>
SoftSudakov::highEnergyRunning(Energy highScale, Energy EWScale,
			       Energy2 s, Energy2 t, Energy2 u,
			       Herwig::EWProcess::Process process,
			       unsigned int iswap) {
  using namespace EWProcess;
  using namespace boost::numeric::ublas;
  using Constants::pi;
  static const Complex I(0,1.0);
  Complex T = getT(s,t), U = getU(s,u);
  matrix<Complex> G1,G2,G3;
  unsigned int numGauge;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    {
      assert(iswap==0);
      numGauge = 4;
      G1 = zero_matrix<Complex>(numGauge,numGauge);
      G2 = zero_matrix<Complex>(numGauge,numGauge);
      G3 = zero_matrix<Complex>(numGauge,numGauge);
      matrix<Complex> gam3 = Gamma3(U,T);
      G3(0,0) += gam3(0,0);
      G3(0,2) += gam3(0,1);
      G3(2,0) += gam3(1,0);
      G3(2,2) += gam3(1,1);
      G3(1,1) += gam3(0,0);
      G3(1,3) += gam3(0,1);
      G3(3,1) += gam3(1,0);
      G3(3,3) += gam3(1,1);
      matrix<Complex> gam2 = Gamma2(U,T);
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
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G2 = Gamma2Singlet();
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,1.0/6.0,T,U);
    break;
  case QQDD:
  case QtQtDD:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G2 = Gamma2Singlet();
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,1.0/6.0,T,U);
    break;
  case QQLL:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3Singlet();
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/2.0,1.0/6.0,T,U);
    break;
  case QQEE:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge); 
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,1.0/6.0,T,U);
    break;
  case UUUU:
  case UUUUiden:
  case tRtRUU:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(2.0/3.0,2.0/3.0,T,U);
    break;
  case UUDD:
  case tRtRDD:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,2.0/3.0,T,U);
    break;
  case UULL:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/2.0,2.0/3.0,T,U);
    break;
  case UUEE:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,2.0/3.0,T,U);
    break;
  case DDDD:
  case DDDDiden:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/3.0,-1.0/3.0,T,U);
    break;
  case DDLL:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/2.0,-1.0/3.0,T,U);
    break;
  case DDEE:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,-1.0/3.0,T,U);
    break;
  case LLLL:
  case LLLLiden:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(-1.0/2.0,-1.0/2.0,T,U);
    break;
  case LLEE:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0,-1.0/2.0,T,U);
    break;
  case EEEE:
  case EEEEiden:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G1(0,0) = Gamma1(-1.0,-1.0,T,U);
    break;
  case QQWW:
    {
      assert(iswap==0);
      numGauge = 5;
      G1 = zero_matrix<Complex>(numGauge,numGauge);
      G3 = zero_matrix<Complex>(numGauge,numGauge);
      Complex gam3s = Gamma3Singlet()(0,0);
      for (unsigned int i=0; i<5; i++) {
	G3(i,i) = gam3s;
	G1(i,i) = Gamma1(1.0/6.0);
      }
      G2 = Gamma2w(U,T);
    }
    break;
  case QQPhiPhi:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = Gamma3Singlet();
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(1.0/2.0,1.0/6.0,T,U);
    break;
  case QQWG:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -17.0/6.0*I*pi + 3.0/2.0*(U+T);
    G2(0,0) = -7.0/4.0*I*pi + (U+T);
    G1(0,0) = Gamma1(1.0/6.0);
    break;
  case QQBG:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G2(0,0) = -3.0/4.0*I*pi;
    G1(0,0) = Gamma1(1.0/6.0);
    break;
  case QQGG:
  case QtQtGG:
    {
      numGauge = 3;
      G1 = zero_matrix<Complex>(numGauge,numGauge);
      G2 = zero_matrix<Complex>(numGauge,numGauge);
      Complex gam2s,gam1;
      if(iswap==0) {
	G3    = Gamma3g(U,T);
	gam2s = Gamma2Singlet()(0,0);
	gam1  = Gamma1(1.0/6.0);
      }
      else if(iswap==1) {
	G3 = Gamma3gST(U,T);
	gam2s = Gamma2SingletST(T)(0,0);
	gam1  = Gamma1ST(1.0/6.0,T);
      }
      else
	assert(false);
      for (unsigned int i=0; i<3; i++) {
	G2(i,i) = gam2s;
	G1(i,i) = gam1;
      }
    }
    break;
  case LLWW:
    assert(iswap==0);
    numGauge = 5;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    for (unsigned int i=0; i<5; i++) {
      G1(i,i) = Gamma1(-1.0/2.0);
    }
    G2 = Gamma2w(U,T);
    break;
  case LLPhiPhi:
    assert(iswap==0);
    numGauge = 2;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = Gamma2(U,T);
    G1(0,0) = G1(1,1) = Gamma1(1.0/2.0,-1.0/2.0,T,U);
    break;
  case UUBB:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(2.0/3.0);
    break;
  case UUPhiPhi:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,2.0/3.0,T,U);
    break;
  case UUBG:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G1(0,0) = Gamma1(2.0/3.0);
    break;
  case DDGG:
  case UUGG:
  case tRtRGG: 
    {
      numGauge = 3;
      G1 = zero_matrix<Complex>(numGauge,numGauge);
      G2 = zero_matrix<Complex>(numGauge,numGauge);
      double Y = process==DDGG ? -1./3. : 2./3.;
      Complex gam1(0.);
      if(iswap==0) {
	G3 = Gamma3g(U,T);
	gam1 = Gamma1(Y);
      }
      else if(iswap==1) {
	G3 = Gamma3gST(U,T);
	gam1 = Gamma1ST(Y,T);
      }
      else
	assert(false);
      for (unsigned int i=0; i<3; i++) {
	G1(i,i) = gam1;
      }
    }
    break;
  case DDBB:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G1(0,0) = Gamma1(-1.0/3.0);
    break;
  case DDPhiPhi:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = Gamma3Singlet()(0,0);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,-1.0/3.0,T,U);
    break;
  case DDBG:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G3(0,0) = -4.0/3.0*I*pi + 3.0/2.0*(U+T-I*pi);
    G1(0,0) = Gamma1(-1.0/3.0);
    break;
  case EEBB:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G1(0,0) = Gamma1(-1.0);
    break;
  case EEPhiPhi:
    assert(iswap==0);
    numGauge = 1;
    G1 = zero_matrix<Complex>(numGauge,numGauge);
    G2 = zero_matrix<Complex>(numGauge,numGauge);
    G3 = zero_matrix<Complex>(numGauge,numGauge);
    G2(0,0) = Gamma2Singlet()(0,0);
    G1(0,0) = Gamma1(1.0/2.0,-1.0,T,U);
    break;
  default:
    assert(false);
  }
  return evaluateSoft(G3,G2,G1,highScale,EWScale,true);
}

unsigned int SoftSudakov::numberGauge(Herwig::EWProcess::Process process) {
  using namespace EWProcess;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    return 4;
  case QQUU:
  case QtQtUU:
  case QQtRtR:
    return 2;
  case QQDD:
  case QtQtDD:
    return 2;
  case QQLL:
    return 2;
  case QQEE:
    return 1;
  case UUUU:
  case UUUUiden:
  case tRtRUU:
    return 2;
  case UUDD:
  case tRtRDD:
    return 2;
  case UULL:
    return 1;
  case UUEE:
    return 1;
  case DDDD:
  case DDDDiden:
    return 2;
  case DDLL:
    return 1;
  case DDEE:
    return 1;
  case LLLL:
  case LLLLiden:
    return 2;
  case LLEE:
    return 1;
  case EEEE:
  case EEEEiden:
    return 1;
  case QQWW:
    return 5;
  case QQPhiPhi:
    return 2;
  case QQWG:
    return 1;
  case QQBG:
    return 1;
  case QQGG:
  case QtQtGG:
    return 3;
  case LLWW:
    return 5;
  case LLPhiPhi:
    return 2;
  case UUBB:
    return 1;
  case UUPhiPhi:
    return 1;
  case UUBG:
    return 1;
  case UUGG:
  case tRtRGG:
    return 3;
  case DDBB:
    return 1;
  case DDPhiPhi:
    return 1;
  case DDBG:
    return 1;
  case DDGG:
    return 3;
  case EEBB:
    return 1;
  case EEPhiPhi:
    return 1;
  default:
    assert(false);
  }
}

unsigned int SoftSudakov::numberBrokenGauge(Herwig::EWProcess::Process process) {
  using namespace EWProcess;
  switch (process) {
  case QQQQ:
  case QQQQiden:
  case QtQtQQ:
    return 12;
  case QQUU:
  case QtQtUU:
  case QQtRtR:
    return 4;
  case QQDD:
  case QtQtDD:
    return 4;
  case QQLL:
    return 6;
  case QQEE:
    return 2;
  case UUUU:
  case UUUUiden:
  case tRtRUU:
    return 2;
  case UUDD:
  case tRtRDD:
    return 2;
  case UULL:
    return 2;
  case UUEE:
    return 1;
  case DDDD:
  case DDDDiden:
    return 2;
  case DDLL:
    return 2;
  case DDEE:
    return 1;
  case LLLL:
  case LLLLiden:
    return 6;
  case EEEE:
  case EEEEiden:
    return 1;
  case QQWW:
    return 20;
  case QQPhiPhi:
    return 14;
  case QQWG:
    return 6;
  case QQBG:
    return 4;
  case QQGG:
  case QtQtGG:
    return 6;
  case LLWW:
    return 20;
  case LLPhiPhi:
    return 14;
  case UUBB:
    return 4;
  case UUPhiPhi:
    return 5;
  case UUBG:
    return 2;
  case UUGG:
  case tRtRGG:
    return 3;
  case DDBB:
    return 4;
  case DDPhiPhi:
    return 5;
  case DDBG:
    return 2;
  case DDGG:
    return 3;
  case EEBB:
    return 4;
  case EEPhiPhi:
    return 5;
  default:
    assert(false);
  }
}
