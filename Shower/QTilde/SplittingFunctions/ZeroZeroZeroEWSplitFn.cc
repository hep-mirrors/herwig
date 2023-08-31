// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZeroZeroZeroEWSplitFn class.
//

#include "ZeroZeroZeroEWSplitFn.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Models/StandardModel/SMFFHVertex.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr ZeroZeroZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr ZeroZeroZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void ZeroZeroZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gw_ << _theSM << _couplingValueIm << _couplingValueRe;
}

void ZeroZeroZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gw_ >> _theSM >> _couplingValueIm >> _couplingValueRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<ZeroZeroZeroEWSplitFn,SplittingFunction>
describeHerwigZeroZeroZeroEWSplitFn("Herwig::ZeroZeroZeroEWSplitFn", "HwShower.so");

void ZeroZeroZeroEWSplitFn::Init() {

  static ClassDocumentation<ZeroZeroZeroEWSplitFn> documentation
    ("The ZeroZeroZeroEWSplitFn class implements the splitting A->A'A''");

  static Parameter<ZeroZeroZeroEWSplitFn, double> interfaceCouplingValueIm
    ("CouplingValue.Im",
     "The numerical value (imaginary part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroZeroEWSplitFn::_couplingValueIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<ZeroZeroZeroEWSplitFn, double> interfaceCouplingValueRe
    ("CouplingValue.Re",
     "The numerical value (real part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroZeroEWSplitFn::_couplingValueRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}

void ZeroZeroZeroEWSplitFn::doinit() {
  SplittingFunction::doinit();
  // set up parameters
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  gw_ = 1./sqrt(sw2);
  // SM cast
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}

void ZeroZeroZeroEWSplitFn::getCouplings(Complex & g, const IdList & ids) const {
  if(ids[0]->iSpin()==PDT::Spin0 && ids[0]->iSpin()==PDT::Spin0 && ids[0]->iSpin()==PDT::Spin0) {
    // BSM cases, where the numerical value of the couplings are
    // expected to be fed into the splitting functions
    if(_couplingValueIm!=0.||_couplingValueRe!=0.) {
      double e  = sqrt(4.*Constants::pi*generator()->standardModel()->
                  alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
      g = Complex(_couplingValueRe,_couplingValueIm)/e;
    }
    // SM case
    else {
      // running masses
      Energy mW = getParticleData(ParticleID::Wplus)->mass();
      Energy mH = getParticleData(ParticleID::h0)->mass();
      g = Complex(0.,1.5*gw_*sqr(mH)/mW/GeV)*10.;
    }
  }
  else
    assert(false);
}

void ZeroZeroZeroEWSplitFn::getCouplings(Complex & g, const IdList & ids,
    const Energy2 t) const {
  if(ids[0]->iSpin()==PDT::Spin0 && ids[0]->iSpin()==PDT::Spin0 && ids[0]->iSpin()==PDT::Spin0) {
    // BSM cases, where the numerical value of the couplings are
    // expected to be fed into the splitting functions
    if(_couplingValueIm!=0.||_couplingValueRe!=0.) {
      double e  = sqrt(4.*Constants::pi*generator()->standardModel()->alphaEM(t));
      g = Complex(_couplingValueRe,_couplingValueIm)/e;
    }
    // SM case
    else {
      // running masses
      Energy mW = _theSM->mass(t,getParticleData(ParticleID::Wplus));
      Energy mH = _theSM->mass(t,getParticleData(ParticleID::h0));
      g = Complex(0.,1.5*gw_*sqr(mH)/mW/GeV)*10.;
    }
  }
  else
    assert(false);
}

double ZeroZeroZeroEWSplitFn::P(const double, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix &) const {
  Complex ghhh(0.,0.);
  getCouplings(ghhh,ids,t);
  double val = norm(ghhh)/(2.*t)*GeV2;
  if(mass)
    return val;
  else
    assert(false);
}

double ZeroZeroZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex ghhh(0.,0.);
  getCouplings(ghhh,ids);
  double val = norm(ghhh);
  return val;
}

double ZeroZeroZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList &, const bool ,
				    const RhoDMatrix & ) const {
  return 1./(2.*t)*GeV2;
}

double ZeroZeroZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex ghhh(0.,0.);
  getCouplings(ghhh,ids);
  double pre = norm(ghhh);
  switch (PDFfactor) {
  case 0: //OverP
    return pre*z;
  case 1: //OverP/z
    return pre*log(z);
  case 2: //OverP/(1-z)
    return -pre*log(-1.+z);
  case 3: //OverP/[z(1-z)]
    return -pre*(log(1.-z)-log(z));
  default:
    throw Exception() << "ZeroZeroZeroEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double ZeroZeroZeroEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex ghhh(0.,0.);
  getCouplings(ghhh,ids);
  double pre = norm(ghhh);
  switch (PDFfactor) {
  case 0:
    return r/pre;
  case 1: //OverP/z
    return exp(r/pre);
  case 2: //OverP/(1-z)
    return 1.+exp(-r/pre);
  case 3: //OverP/[z(1-z)]
    return exp(r/pre)/(1.+exp(r/pre));
  default:
    throw Exception() << "ZeroZeroZeroEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool ZeroZeroZeroEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3)
    return false;
  if(ids[0]->mass()==0.*GeV)
    return false;
  if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge())
    return false;
  if(ids[0]->iSpin()==PDT::Spin0 && ids[1]->iSpin()==PDT::Spin0 && ids[2]->iSpin()==PDT::Spin0)
    return true;
  else
    return false;
}

vector<pair<int, Complex> >
ZeroZeroZeroEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // scalar so no dependence
  return {{ {0, 1.} }};
}

vector<pair<int, Complex> >
ZeroZeroZeroEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // scalar so no dependence
  assert(false);
  return {{ {0, 1.} }};
}

DecayMEPtr ZeroZeroZeroEWSplitFn::matrixElement(const double, const Energy2 t,
                                             const IdList & ids, const double,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  Complex ghhh(0.,0.);
  getCouplings(ghhh,ids,t);
  double sqrtt = sqrt(t)/GeV;
  (*kernal)(0,0,0) = ghhh/sqrt(2)/sqrtt;
  // return the answer
  return kernal;
}
