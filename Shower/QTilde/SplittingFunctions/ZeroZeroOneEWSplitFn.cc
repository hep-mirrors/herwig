// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZeroZeroOneEWSplitFn class.
//

#include "ZeroZeroOneEWSplitFn.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Models/StandardModel/SMFFHVertex.h"

using namespace Herwig;

IBPtr ZeroZeroOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr ZeroZeroOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void ZeroZeroOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gHHZ_ << _theSM;
}

void ZeroZeroOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gHHZ_ >> _theSM;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<ZeroZeroOneEWSplitFn,SplittingFunction>
describeHerwigZeroZeroOneEWSplitFn("Herwig::ZeroZeroOneEWSplitFn", "HwShower.so");


void ZeroZeroOneEWSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneEWSplitFn> documentation
    ("The ZeroZeroOneEWSplitFn class implements the splitting H->HZ");

}


void ZeroZeroOneEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // coupling
  gHHZ_ = 0.5/(sqrt(sw2)*sqrt(1.-sw2));
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}


void ZeroZeroOneEWSplitFn::getCouplings(double & g, const IdList & ids) const {
  if(ids[0]->id()==ParticleID::h0 && ids[1]->id()==ParticleID::h0 && ids[2]->id()==ParticleID::Z0){
    g = gHHZ_;
  }
  else
    assert(false);
}


double ZeroZeroOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gHHZ(0.);
  getCouplings(gHHZ,ids);
  // massless limit
  double val = 4.*z/(1.-z);
  // massive limits
  if(mass) {
    // running masses of the bosons
    Energy mH_r = _theSM->mass(t,getParticleData(ParticleID::h0));
    Energy mZ_r = _theSM->mass(t,getParticleData(ParticleID::Z0));
    double mHt2 = sqr(mH_r)/t;
    double mZt2 = sqr(mZ_r)/t;
    val += 2.*mZt2*(1.+sqr(z))/sqr(1.-z);
    val -= 4.*mHt2;
  }
  return sqr(gHHZ)*val;
}


double ZeroZeroOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gHHZ(0.);
  getCouplings(gHHZ,ids);
  return sqr(gHHZ)*4./(1.-z);
}


double ZeroZeroOneEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  // massless limit
  double val = z;
  // massive limits
  if(mass) {
    // running masses of the bosons
    Energy mH_r = _theSM->mass(t,getParticleData(ParticleID::h0));
    Energy mZ_r = _theSM->mass(t,getParticleData(ParticleID::Z0));
    double mHt2 = sqr(mH_r)/t;
    double mZt2 = sqr(mZ_r)/t;
    val += 0.5*mZt2*(1.+sqr(z))/(1.-z);
    val -= mHt2*(1.-z);
  }
  return val;
}


double ZeroZeroOneEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gHHZ(0.);
  getCouplings(gHHZ,ids);
  double pre = sqr(gHHZ);
  switch (PDFfactor) {
  case 0:
    return -4.*pre*log(1.-z);
  case 1:
    //return -2.*pre*(1./z+log(1.-z)-log(z));
  case 2:
    //return 2.*pre*(2.*log(z)+(2.*z-1.)/(z*(1.-z))-2.*log(1.-z));
  case 3:
    //return 2.*pre*(1./(1.-z)-1./z-2.*log(1.-z)+2.*log(z));
  default:
    throw Exception() << "ZeroZeroOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double ZeroZeroOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  double gHHZ(0.);
  getCouplings(gHHZ,ids);
  double pre = 4.*sqr(gHHZ);
  switch (PDFfactor) {
  case 0:
    return 1.-exp(-r/pre);
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "ZeroZeroOneEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool ZeroZeroOneEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3)
    return false;
  if(ids[0]->id()==ParticleID::h0 && ids[1]->id()==ParticleID::h0 && ids[2]->id()==ParticleID::Z0)
    return true;
  else
    return false;
}


vector<pair<int, Complex> >
ZeroZeroOneEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


vector<pair<int, Complex> >
ZeroZeroOneEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


DecayMEPtr ZeroZeroOneEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  double gHHZ(0.);
  getCouplings(gHHZ,ids);
  // defining dummies
  Energy mH_r = _theSM->mass(t,getParticleData(ParticleID::h0));
  Energy mZ_r = _theSM->mass(t,getParticleData(ParticleID::Z0));
  double m0t = mH_r/sqrt(t);
  double m1t = mH_r/sqrt(t);
  double m2t = mZ_r/sqrt(t);
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  double sqrtmass = sqrt(sqr(m0t)-sqr(m1t)/z-sqr(m2t)/(1.-z)+1.);
  // assign kernel
  (*kernal)(0,0,0) = 2.*gHHZ*phase*sqrt(z/(1.-z))*sqrtmass;
  (*kernal)(0,0,1) = sqrt(2.)*gHHZ*m2t*((1.+z)/(1.-z)); //2>4
  (*kernal)(0,0,2) = -2.*gHHZ*cphase*sqrt(z/(1.-z))*sqrtmass;
  // return the answer
  return kernal;
}
