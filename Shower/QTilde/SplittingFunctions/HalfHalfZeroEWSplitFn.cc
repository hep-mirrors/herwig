// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfZeroEWSplitFn class.
//

#include "HalfHalfZeroEWSplitFn.h"
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

IBPtr HalfHalfZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HalfHalfZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void HalfHalfZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << ghqq_ << _theSM;
}

void HalfHalfZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> ghqq_ >> _theSM;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<HalfHalfZeroEWSplitFn,SplittingFunction>
describeHerwigHalfHalfZeroEWSplitFn("Herwig::HalfHalfZeroEWSplitFn", "HwShower.so");

void HalfHalfZeroEWSplitFn::Init() {

  static ClassDocumentation<HalfHalfZeroEWSplitFn> documentation
    ("The HalfHalfZeroEWSplitFn class implements the splitting q->qH");

}

void HalfHalfZeroEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  ghqq_ = 1./sqrt(4.*sw2);

  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}

void HalfHalfZeroEWSplitFn::getCouplings(double & gH, const IdList & ids) const {
  if(abs(ids[2]->id())==ParticleID::h0) {
    //get quark masses
    Energy mq;
    if(abs(ids[0]->id())==ParticleID::c)
      mq = getParticleData(ParticleID::c)->mass();
    else if(abs(ids[0]->id())==ParticleID::b)
      mq = getParticleData(ParticleID::b)->mass();
    else if(abs(ids[0]->id())==ParticleID::t)
      mq = getParticleData(ParticleID::t)->mass();
    Energy mW = getParticleData(ParticleID::Wplus)->mass();
    gH = ghqq_*(mq/mW);
  }
  else
    assert(false);
}

void HalfHalfZeroEWSplitFn::getCouplings(double & gH, const IdList & ids, const Energy2 t) const {
  if(abs(ids[2]->id())==ParticleID::h0) {
    //get quark masses
    Energy mq;
    if(abs(ids[0]->id())==ParticleID::c)
      mq = _theSM->mass(t,getParticleData(ParticleID::c));
    else if(abs(ids[0]->id())==ParticleID::b)
      mq = _theSM->mass(t,getParticleData(ParticleID::b));
    else if(abs(ids[0]->id())==ParticleID::t)
      mq = _theSM->mass(t,getParticleData(ParticleID::t));
    Energy mW = getParticleData(ParticleID::Wplus)->mass();
    //Energy mW = _theSM->mass(t,getParticleData(ParticleID::Wplus));
    gH = ghqq_*(mq/mW);
  }
  else
    assert(false);
}

double HalfHalfZeroEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gH(0.);
  getCouplings(gH,ids,t);
  double val = (1.-z);
  Energy mq, mH;
  //get masses
  if(mass) {
    mq = ids[0]->mass();
    mH = ids[2]->mass();
  }
  else { // to assure the particle mass in non-zero
    if(abs(ids[0]->id())==ParticleID::c)
      mq = getParticleData(ParticleID::c)->mass();
    else if(abs(ids[0]->id())==ParticleID::b)
      mq = getParticleData(ParticleID::b)->mass();
    else if(abs(ids[0]->id())==ParticleID::t)
      mq = getParticleData(ParticleID::t)->mass();
    mH = getParticleData(ParticleID::h0)->mass();
  }
  val += (4.*sqr(mq) - sqr(mH))/(t*(1. - z)*z);
  val *= sqr(gH);
  return colourFactor(ids)*val;
}


double HalfHalfZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gH(0.);
  getCouplings(gH,ids);
  return sqr(gH)*colourFactor(ids)*(1.-z);
}

double HalfHalfZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double gH(0.);
  getCouplings(gH,ids,t);
  double val = 1.;
  Energy mq, mH;
  if(mass) {
    mq = ids[0]->mass();
    mH = ids[2]->mass();
  }
  else { // to assure the particle mass in non-zero
    if(abs(ids[0]->id())==ParticleID::c)
      mq = getParticleData(ParticleID::c)->mass();
    else if(abs(ids[0]->id())==ParticleID::b)
      mq = getParticleData(ParticleID::b)->mass();
    else if(abs(ids[0]->id())==ParticleID::t)
      mq = getParticleData(ParticleID::t)->mass();
    mH = getParticleData(ParticleID::h0)->mass();
  }
  val += (4.*sqr(mq) - sqr(mH))/(t*(1. - z)*z);
  return val;
}

double HalfHalfZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gH(0.);
  getCouplings(gH,ids);
  double pre = colourFactor(ids)*sqr(gH);
  switch (PDFfactor) {
  case 0: //OverP
    return pre*(z-sqr(z)/2.);
  case 1: //OverP/z
    return  pre*(log(z)-z);
  case 2: //OverP/(1-z)
    return  pre*z;
  case 3: //OverP/[z(1-z)]
    return  pre*log(z);
  default:
    throw Exception() << "HalfHalfZeroEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfHalfZeroEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  double gH(0.);
  getCouplings(gH,ids);
  double pre = colourFactor(ids)*sqr(gH);
  switch (PDFfactor) {
  case 0:
    return 1. - sqrt(1. - 2.*r/pre);
  case 1: //OverP/z
  case 2: //OverP/(1-z)
    return  r/pre;
  case 3: //OverP/[z(1-z)]
    return  exp(r/pre);
  default:
    throw Exception() << "HalfHalfZeroEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfHalfZeroEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[2]->id()==ParticleID::h0) {
    if(ids[0]->id()==ids[1]->id() && (ids[0]->id()==4 || ids[0]->id()==5 || ids[0]->id()==6))
      return true;
  }
  return false;
}

vector<pair<int, Complex> >
HalfHalfZeroEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

vector<pair<int, Complex> >
HalfHalfZeroEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

DecayMEPtr HalfHalfZeroEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  //get masses
  Energy mq, mH;
  if(abs(ids[0]->id())==ParticleID::c)
    mq = getParticleData(ParticleID::c)->mass();
  else if(abs(ids[0]->id())==ParticleID::b)
    mq = getParticleData(ParticleID::b)->mass();
  else if(abs(ids[0]->id())==ParticleID::t)
    mq = getParticleData(ParticleID::t)->mass();
  mH = getParticleData(ParticleID::h0)->mass();
  double gH(0.);
  getCouplings(gH,ids,t);
  double mqt = mq/sqrt(t);
  double mHt = mH/sqrt(t);
  double num1 = gH*(1.+z)*mqt;
  double num2 = gH*sqrt(-sqr(mqt)*(1.-z) - sqr(mHt)*z + z*(1.-z)*(sqr(mqt)+z*(1.-z))); //watch this
  double dnum = sqrt(2.)*sqrt((1.-z)*sqr(z));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  (*kernal)(0,0,0) = num1/dnum;
  (*kernal)(0,1,0) = cphase*num2/dnum;
  (*kernal)(1,0,0) = -phase*num2/dnum;
  (*kernal)(1,1,0) = num1/dnum;
  // return the answer
  return kernal;
}
