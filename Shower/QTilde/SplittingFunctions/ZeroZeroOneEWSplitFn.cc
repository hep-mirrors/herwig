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
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr ZeroZeroOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr ZeroZeroOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void ZeroZeroOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << _couplingValueIm << _couplingValueRe;
}

void ZeroZeroOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> _couplingValueIm >> _couplingValueRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<ZeroZeroOneEWSplitFn,SplittingFunction>
describeHerwigZeroZeroOneEWSplitFn("Herwig::ZeroZeroOneEWSplitFn", "HwShower.so");


void ZeroZeroOneEWSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneEWSplitFn> documentation
    ("The ZeroZeroOneEWSplitFn class implements purly beyond SM electroweak splittings H->H'V");

  static Parameter<ZeroZeroOneEWSplitFn, double> interfaceCouplingValueIm
    ("CouplingValue.Im",
     "The numerical value (imaginary part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroOneEWSplitFn::_couplingValueIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<ZeroZeroOneEWSplitFn, double> interfaceCouplingValueRe
    ("CouplingValue.Re",
     "The numerical value (real part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroOneEWSplitFn::_couplingValueRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}


void ZeroZeroOneEWSplitFn::doinit() {
  SplittingFunction::doinit();
}


void ZeroZeroOneEWSplitFn::getCouplings(Complex & g, const IdList &) const {
  g = Complex(_couplingValueRe,_couplingValueIm);
}


double ZeroZeroOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  Complex ghhv(0.,0.);
  getCouplings(ghhv,ids);
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  // the splitting in the massless limit
  double val = z*(rho00+rho22)/(1.-z);
  // the massive limit
  if(mass){
    // get the running mass
    double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
    double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
    double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
    val += (m2t2*sqr(1.+z)*rho11)/(2.*sqr(1.-z))+((-(m1t2*(1.-z))
        -m2t2*z+m0t2*(1.-z)*z)*(rho00+rho22))/sqr(1.-z);
  }
  return norm(ghhv)*val;
}


double ZeroZeroOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex ghhv(0.,0.);
  getCouplings(ghhv,ids);
  return norm(ghhv)*z/(1.-z);
}


double ZeroZeroOneEWSplitFn::ratioP(const double z, const Energy2 t,
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
    val += ((1.-z)*(m2t2*sqr(1.+z)*rho11+2*(m1t2*(-1.+z)
        -(m2t2+m0t2*(-1.+z))*z)*(rho00+rho22)))/(2.*sqr(-1.+z)*z);
  }
  return val;
}


double ZeroZeroOneEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex ghhv(0.,0.);
  getCouplings(ghhv,ids);
  double pre = norm(ghhv);
  switch (PDFfactor) {
  case 0:
    return -pre*(z+log(1.-z));
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "ZeroZeroOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}


double ZeroZeroOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex ghhv(0.,0.);
  getCouplings(ghhv,ids);
  double pre = norm(ghhv);
  switch (PDFfactor) {
  case 0:
    return 1.-exp(-(1+r/pre));;
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
  if(_couplingValueIm==0.&&_couplingValueRe==0.)
    return false;
  if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge())
    return false;
  if(ids[0]->iSpin()==PDT::Spin0 && ids[1]->iSpin()==PDT::Spin1 && ids[2]->iSpin()==PDT::Spin0)
    return true;
  else if(ids[0]->iSpin()==PDT::Spin0 && ids[1]->iSpin()==PDT::Spin0 && ids[2]->iSpin()==PDT::Spin1)
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
  Complex ghhv(0.,0.);
  getCouplings(ghhv,ids);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
  double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
  double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  double sqrtmass = sqrt(m0t2*z*(1.-z)-m1t2*(1.-z)-m2t2*z+z*(1.-z));
  // assign kernel
  (*kernal)(0,0,0) = -phase*ghhv*sqrtmass/(1.-z);        // 111
  (*kernal)(1,0,0) = -sqrt(m2t2)*(1.+z)/sqrt(2.*(1.-z)); // 211 -> 411
  (*kernal)(2,0,0) = cphase*ghhv*sqrtmass/(1.-z);        // 311
  // return the answer
  return kernal;
}
