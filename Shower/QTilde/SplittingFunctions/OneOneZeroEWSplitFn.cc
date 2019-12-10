// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneZeroEWSplitFn class.
//

#include "OneOneZeroEWSplitFn.h"
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

IBPtr OneOneZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneOneZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneOneZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gWWH_ << gZZH_ << _theSM;
}

void OneOneZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gWWH_ >> gZZH_ >> _theSM;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneOneZeroEWSplitFn,SplittingFunction>
describeHerwigOneOneZeroEWSplitFn("Herwig::OneOneZeroEWSplitFn", "HwShower.so");


void OneOneZeroEWSplitFn::Init() {

  static ClassDocumentation<OneOneZeroEWSplitFn> documentation
    ("The OneOneZeroEWSplitFn class implements the splittings W->WH and Z->ZH");

}


void OneOneZeroEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // W -> W H coupling      g = e/sin theta_W
  gWWH_ = 1./sqrt(sw2);
  // Z -> Z H coupling
  gZZH_ = 1./(sqrt(sw2)*sqrt(1.-sw2));
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}


void OneOneZeroEWSplitFn::getCouplings(double & g, const IdList & ids) const {
  if(abs(ids[0]->id())==ParticleID::Wplus){
    g = gWWH_;
  }
  else if(ids[0]->id()==ParticleID::Z0){
    g = gZZH_;
  }
  else
    assert(false);
}


double OneOneZeroEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gvvh(0.);
  getCouplings(gvvh,ids);
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  // the splitting in the massless limit
  double m0t = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
    //val += (2./sqr(z))*sqr(m0t)*(rho11+sqr(z)*(rho00+rho22));
  double val = 2.*sqr(m0t)*(rho00+rho22);
  return sqr(gvvh)*val;
}


double OneOneZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gvvh(0.);
  getCouplings(gvvh,ids);
  return sqr(gvvh)*2.;
}


double OneOneZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  double m0t = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
  double val = sqr(m0t)*(rho00+rho22);
  return val;
}


double OneOneZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gvvh(0.);
  getCouplings(gvvh,ids);
  double pre = sqr(gvvh);
  switch (PDFfactor) {
  case 0:
    return pre*2.*z;
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneZeroEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneOneZeroEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  double gvvh(0.);
  getCouplings(gvvh,ids);
  double pre = sqr(gvvh);
  switch (PDFfactor) {
  case 0:
    return r/(pre*2.);
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneZeroEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneZeroEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3)
    return false;
  if(ids[0]->id()!=ids[1]->id())
    return false;
  if(abs(ids[0]->id())==ParticleID::Wplus && ids[2]->id()==ParticleID::h0)
    return true;
  else if(ids[0]->id()==ParticleID::Z0 && ids[2]->id()==ParticleID::h0)
    return true;
  else
    return false;
}


vector<pair<int, Complex> >
OneOneZeroEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


vector<pair<int, Complex> >
OneOneZeroEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


DecayMEPtr OneOneZeroEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  double gvvh(0.);
  getCouplings(gvvh,ids);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  double r2 = sqrt(2.);
  double m0t = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
  double m2t = _theSM->mass(t,getParticleData(ids[2]->id()))/sqrt(t);
  double sqrtmass = sqrt(sqr(m0t)-sqr(m0t)/z-sqr(m2t)/(1.-z)+1.);
  // assign kernel
  (*kernal)(0,0,0) = -gvvh*r2*m0t; // 111
  (*kernal)(0,1,0) = 0.; // 121
  (*kernal)(0,2,0) = 0.; // 131
  (*kernal)(1,0,0) = 0.; // 211
  (*kernal)(1,1,0) = 0.; // 221 > 441
  (*kernal)(1,2,0) = 0.; // 231
  (*kernal)(2,0,0) = 0.; // 311
  (*kernal)(2,1,0) = 0.; // 321 > 341
  (*kernal)(2,2,0) = -gvvh*r2*m0t; // 331
  // return the answer
  return kernal;
}
