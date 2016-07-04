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
