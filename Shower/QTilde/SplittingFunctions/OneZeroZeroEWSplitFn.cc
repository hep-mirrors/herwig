// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneZeroZeroEWSplitFn class.
//

#include "OneZeroZeroEWSplitFn.h"
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

IBPtr OneZeroZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneZeroZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneZeroZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << _couplingValue;
}

void OneZeroZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> _couplingValue;;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneZeroZeroEWSplitFn,SplittingFunction>
describeHerwigOneZeroZeroEWSplitFn("Herwig::OneZeroZeroEWSplitFn", "HwShower.so");


void OneZeroZeroEWSplitFn::Init() {

  static ClassDocumentation<OneZeroZeroEWSplitFn> documentation
    ("The OneZeroZeroEWSplitFn class implements purly beyond SM electroweak splittings V->HH'");

  static Parameter<OneZeroZeroEWSplitFn, double> interfaceCouplingValue
    ("CouplingValue",
     "The numerical value of the splitting coupling to be imported for BSM splittings",
     &OneZeroZeroEWSplitFn::_couplingValue, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}


void OneZeroZeroEWSplitFn::doinit() {
  SplittingFunction::doinit();
}


void OneZeroZeroEWSplitFn::getCouplings(double & g, const IdList &) const {
  g = _couplingValue;
}


double OneZeroZeroEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gvhh(0.);
  getCouplings(gvhh,ids);
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  // the splitting in the massless limit
  double val = (1.-z)*z*(rho00+rho22);
  // the massive limit
  if(mass){
    // get the running mass
    double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
    double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
    double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
    val += (m0t2*sqr(-1.+2.*z)*rho11)/2.+(-(m1t2*(1.-z))-m2t2*z+m0t2*(1.-z)*z)
        *(rho00+rho22);
  }
  return sqr(gvhh)*val;
}


double OneZeroZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gvhh(0.);
  getCouplings(gvhh,ids);
  return sqr(gvhh)*(1.-z)*z;
}


double OneZeroZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  // ratio in the massless limit
  double val = rho00+rho22;
  // the massive limit
  if(mass){
    // get the running mass
    double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
    double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
    double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
    val += ((m0t2*sqr(-1.+2.*z)*rho11)/2.+(-(m1t2*(1.-z))-m2t2*z+m0t2*(1.-z)*z)
        *(rho00+rho22))/((1.-z)*z);
  }
  return val;
}


double OneZeroZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gvhh(0.);
  getCouplings(gvhh,ids);
  double pre = sqr(gvhh);
  switch (PDFfactor) {
  case 0:
    return pre/6.*(3.-2.*z)*sqr(z);
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneZeroZeroEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}


double OneZeroZeroEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  double gvhh(0.);
  getCouplings(gvhh,ids);
  double pre = sqr(gvhh);
  switch (PDFfactor) {
  case 0:
    return (1.-pre/pow(-pow(pre,3)+12.*pow(pre,2)*r+2.*sqrt(6.)
            *sqrt(-(pow(pre,4)*(pre-6.*r)*r)),1./3.)- pow(-pow(pre,3)+12.
            *pow(pre,2)*r+2.*sqrt(6)*sqrt(-(pow(pre,4)*(pre-6.*r)*r)),1./3.)/pre)/2.;
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneZeroZeroEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}


bool OneZeroZeroEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3)
    return false;
  if(_couplingValue==0.)
    return false;
  if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge())
    return false;
  bool isGVB = (abs(ids[0]->id())==ParticleID::Wplus
             || abs(ids[0]->id())==ParticleID::Z0
             || abs(ids[0]->id())==ParticleID::gamma);
  if(isGVB && ids[1]->iSpin()==PDT::Spin0 && ids[2]->iSpin()==PDT::Spin0)
    return true;
  else
    return false;
}


vector<pair<int, Complex> >
OneZeroZeroEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


vector<pair<int, Complex> >
OneZeroZeroEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


DecayMEPtr OneZeroZeroEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  double gvhh(0.);
  getCouplings(gvhh,ids);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
  double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
  double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  double sqrtmass = sqrt(m0t2*z*(1.-z)-m1t2*(1.-z)-m2t2*z+z*(1.-z));
  // assign kernel
  (*kernal)(0,0,0) = -cphase*sqrtmass;                        // 111
  (*kernal)(1,0,0) =  sqrt(m0t2)*(1.-2.*z)/sqrt(2.);         // 211 -> 411
  (*kernal)(2,0,0) =  phase*sqrtmass;                         // 311
  // return the answer
  return kernal;
}
