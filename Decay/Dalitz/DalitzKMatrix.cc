// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzKMatrix class.
//

#include "DalitzKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzKMatrix::persistentOutput(PersistentOStream & os) const {
  os << kMatrix_ << channel_ << imat_ << ounit(sc_,GeV2) << expType_ << beta_ << coeffs_;
}

void DalitzKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> kMatrix_ >> channel_ >> imat_ >> iunit(sc_,GeV2) >> expType_ >> beta_ >> coeffs_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzKMatrix,DalitzResonance>
describeHerwigDalitzKMatrix("Herwig::DalitzKMatrix", "HwDalitzDecay.so");

void DalitzKMatrix::Init() {

  static ClassDocumentation<DalitzKMatrix> documentation
    ("The DalitzKMatrix class allows the use of \f$K\f$-matrices in Dalitz decays");

}

Complex DalitzKMatrix::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  Energy2 s = sqr(mAB);
  double sHat = (s-sc_)/GeV2;
  // construct the p-vector
  ublas::vector<Complex> pVector(coeffs_.size());
  // compute the terms
  for(unsigned int ix=0;ix<coeffs_.size();++ix) {
    Complex val(0.);
    // first the pole piece
    for(unsigned int iy=0;iy<kMatrix_->poles().size();++iy) {
      Complex piece = GeV*beta_[iy]*kMatrix_->poleCouplings()[iy][ix]/kMatrix_->poles()[iy];
      for(unsigned int iz=0;iz<kMatrix_->poles().size();++iz) {
	if(iz==iy) continue;
	piece *= 1. - s/kMatrix_->poles()[iz];
      }
      val +=piece;
    }
    Complex fact=exp(Complex(0.,coeffs_[ix].first));
    for(unsigned int iz=0;iz<kMatrix_->poles().size();++iz)
      fact *= 1. - s/kMatrix_->poles()[iz];
    // then the polynomial piece
    double poly=coeffs_[ix].second[0];
    if(expType_==0) {
      for(unsigned int iz=1;iz<coeffs_[ix].second.size();++iz)
	poly += coeffs_[ix].second[iz]*pow(sHat,iz);
    }
    else {
      poly *= (GeV2-sc_)/(s-sc_);
    }
    // store the answer
    pVector[ix]=val+fact*poly;
  }
  ublas::vector<Complex> amps = kMatrix_->amplitudes(s,pVector,true);
  return amps[channel_];
}

void DalitzKMatrix::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << imat_ << " " << channel_ << " " << sc_/GeV2;
  for(unsigned int ix=0;ix<beta_.size();++ix)
    output << " " << abs(beta_[ix]) << " " << arg(beta_[ix]);
  for(unsigned int ix=0;ix<coeffs_.size();++ix) {
    output << " " << coeffs_[ix].second.size();
    for(unsigned int iy=0;iy<coeffs_[ix].second.size();++iy)
      output << " " << coeffs_[ix].second[iy];
    output << " " << coeffs_[ix].first;
  }
}
