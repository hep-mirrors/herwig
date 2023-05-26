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
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr HalfHalfZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HalfHalfZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void HalfHalfZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << ghqq_ << _theSM << _couplingValue0Im << _couplingValue0Re << _couplingValue1Im << _couplingValue1Re;
}

void HalfHalfZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> ghqq_ >> _theSM >> _couplingValue0Im >> _couplingValue0Re >> _couplingValue1Im >> _couplingValue1Re;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<HalfHalfZeroEWSplitFn,SplittingFunction>
describeHerwigHalfHalfZeroEWSplitFn("Herwig::HalfHalfZeroEWSplitFn", "HwShower.so");

void HalfHalfZeroEWSplitFn::Init() {

  static ClassDocumentation<HalfHalfZeroEWSplitFn> documentation
    ("The HalfHalfZeroEWSplitFn class implements the splitting q->qH");

  static Parameter<HalfHalfZeroEWSplitFn, double> interfaceCouplingValueCP0Im
    ("CouplingValue.CP0.Im",
     "The numerical value (imaginary part) of the CP-even splitting coupling to be imported for BSM splittings",
     &HalfHalfZeroEWSplitFn::_couplingValue0Im, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfZeroEWSplitFn, double> interfaceCouplingValueCP0Re
    ("CouplingValue.CP0.Re",
     "The numerical value (real part) of the CP-even splitting coupling to be imported for BSM splittings",
     &HalfHalfZeroEWSplitFn::_couplingValue0Re, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfZeroEWSplitFn, double> interfaceCouplingValueCP1Im
    ("CouplingValue.CP1.Im",
     "The numerical value (imaginary part) of the CP-odd splitting coupling to be imported for BSM splittings",
     &HalfHalfZeroEWSplitFn::_couplingValue1Im, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfZeroEWSplitFn, double> interfaceCouplingValueCP1Re
    ("CouplingValue.CP1.Re",
     "The numerical value (real part) of the CP-odd splitting coupling to be imported for BSM splittings",
     &HalfHalfZeroEWSplitFn::_couplingValue1Re, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}

void HalfHalfZeroEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  ghqq_ = 1./sqrt(4.*sw2);

  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}

void HalfHalfZeroEWSplitFn::getCouplings(Complex& gH0, Complex& gH1, const IdList & ids) const {
  if(ids[2]->iSpin()==PDT::Spin0) {
    if(_couplingValue0Re!=0||_couplingValue0Im!=0||_couplingValue1Re!=0||_couplingValue1Im!=0) {
      double e = sqrt(4.*Constants::pi*generator()->standardModel()
                ->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
      // a factor e is factored out since its already accounted for
      gH0 = Complex(_couplingValue0Re,_couplingValue0Im)/e;
      gH1 = Complex(_couplingValue1Re,_couplingValue1Im)/e;
    }
    else {
      //get quark masses
      Energy m0 = getParticleData(abs(ids[0]->id()))->mass();
      Energy mW = getParticleData(ParticleID::Wplus)->mass();
      gH0 = ghqq_*(m0/mW);
    }
  }
  else
    assert(false);
}

void HalfHalfZeroEWSplitFn::getCouplings(Complex& gH0, Complex& gH1, const IdList & ids, const Energy2 t) const {
  if(ids[2]->iSpin()==PDT::Spin0) {
    if(_couplingValue0Re!=0||_couplingValue0Im!=0||_couplingValue1Re!=0||_couplingValue1Im!=0) {
      Energy mZ = _theSM->mass(t,getParticleData(ParticleID::Z0));
      double e  = sqrt(4.*Constants::pi*generator()->standardModel()->alphaEM(sqr(mZ)));
      gH0 = Complex(_couplingValue0Re,_couplingValue0Im)/e;
      gH1 = Complex(_couplingValue1Re,_couplingValue1Im)/e;
    }
    else {
      //get quark masses
      Energy m0 = _theSM->mass(t,getParticleData(abs(ids[0]->id())));
      Energy mW = getParticleData(ParticleID::Wplus)->mass();
      //Energy mW = _theSM->mass(t,getParticleData(ParticleID::Wplus));
      gH0 = ghqq_*(m0/mW);
    }
  }
  else
    assert(false);
}

double HalfHalfZeroEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids,t);
  double val = (abs(rho(1,1))*norm(gH0+gH1)+abs(rho(0,0))*norm(gH0-gH1))*(1-z);
  Energy m0, m1, m2;
  //get masses
  if(mass) {
    m0 = ids[0]->mass();
    m1 = ids[1]->mass();
    m2  = ids[2]->mass();
  }
  else { // to assure the particle mass in non-zero
    m0 = getParticleData(abs(ids[0]->id()))->mass();
    m1 = getParticleData(abs(ids[1]->id()))->mass();
    m2  = getParticleData(abs(ids[2]->id()))->mass();
  }
  double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);
  val += -(abs(rho(1,1))*norm(gH0+gH1)+abs(rho(0,0))*norm(gH0-gH1))*sqr(m2t)
         +(abs(rho(1,1))+abs(rho(0,0)))*(norm(gH0)*sqr(m0t+m1t)+norm(gH1)*sqr(m0t-m1t))
         +2*(abs(rho(1,1))-abs(rho(0,0)))*(gH0*conj(gH1)).real()*((1-2*z)*sqr(m0t)+sqr(m1t));
  return colourFactor(ids)*val;
}


double HalfHalfZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids);
  return colourFactor(ids)*(norm(gH0)+norm(gH1)+2*abs((gH0*gH1).real())+2*abs((gH0*gH1).imag()))*(1.-z);
}

double HalfHalfZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids,t);
  double val = (abs(rho(1,1))*norm(gH0+gH1)+abs(rho(0,0))*norm(gH0-gH1))*(1-z);
  Energy m0, m1, m2;
  if(mass) {
    m0 = ids[0]->mass();
    m1 = ids[1]->mass();
    m2  = ids[2]->mass();
  }
  else { // to enssure the particle mass in non-zero
    m0 = getParticleData(abs(ids[0]->id()))->mass();
    m1 = getParticleData(abs(ids[1]->id()))->mass();
    m2  = getParticleData(ids[2]->id())->mass();
  }
  double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);
  val += -(abs(rho(1,1))*norm(gH0+gH1)+abs(rho(0,0))*norm(gH0-gH1))*sqr(m2t)
         +(abs(rho(1,1))+abs(rho(0,0)))*(norm(gH0)*sqr(m0t+m1t)+norm(gH1)*sqr(m0t-m1t))
         +2*(abs(rho(1,1))-abs(rho(0,0)))*(gH0*conj(gH1)).real()*((1-2*z)*sqr(m0t)+sqr(m1t));
  val /= (1.-z)*(norm(gH0)+norm(gH1)+2*abs((gH0*gH1).real())+2*abs((gH0*gH1).imag()));
  return val;
}

double HalfHalfZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids);
  double pre = colourFactor(ids)*(norm(gH0)+norm(gH1)+2*abs((gH0*gH1).real())+2*abs((gH0*gH1).imag()));
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
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids);
  double pre = colourFactor(ids)*(norm(gH0)+norm(gH1)+2*abs((gH0*gH1).real())+2*abs((gH0*gH1).imag()));
  switch (PDFfactor) {
  case 0:
    return 1.-sqrt(1.-2.*r/pre);
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
  if(ids[2]->iSpin()==PDT::Spin0 && _couplingValue0Re==0 && _couplingValue0Im==0 && _couplingValue1Re==0 && _couplingValue1Im==0) {
    if(ids[0]->id()==ids[1]->id() && (ids[0]->id()==4 || ids[0]->id()==5 || ids[0]->id()==6)) return true;
  }
  else if(ids[2]->iSpin()==PDT::Spin0 && !(_couplingValue0Re==0 && _couplingValue0Im==0 && _couplingValue1Re==0 && _couplingValue1Im==0)) {
    if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge()) return false;
    if((abs(ids[0]->id())>=1 && abs(ids[0]->id())<=6) && (abs(ids[1]->id())>=1 && abs(ids[1]->id())<=6)) return true;
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
  Energy m0 = getParticleData(abs(ids[0]->id()))->mass();
  Energy m1 = getParticleData(abs(ids[1]->id()))->mass();
  Energy m2  = getParticleData(abs(ids[2]->id()))->mass();
  Complex gH0(0.,0.);
  Complex gH1(0.,0.);
  getCouplings(gH0,gH1,ids,t);
  Energy den = sqrt(2*t);
  Energy pt  = sqrt(z*(1-z)*t+z*(1-z)*sqr(m0)-(1-z)*sqr(m1)-z*sqr(m2));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  (*kernal)(0,0,0) = (gH0*(z*m0+m1)+gH1*(z*m0-m1))/sqrt(z)/den;
  (*kernal)(0,1,0) = cphase*(gH0-gH1)*pt/sqrt(z)/den;
  (*kernal)(1,0,0) = -phase*(gH0+gH1)*pt/sqrt(z)/den;
  (*kernal)(1,1,0) = (gH0*(z*m0+m1)-gH1*(z*m0-m1))/sqrt(z)/den;
  // return the answer
  return kernal;
}
