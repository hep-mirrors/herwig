// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneEWSplitFn class.
//

#include "HalfHalfOneEWSplitFn.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr HalfHalfOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HalfHalfOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void HalfHalfOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gZ_ << gWL_ << _couplingValueLeftIm << _couplingValueLeftRe << _couplingValueRightIm << _couplingValueRightRe;
}

void HalfHalfOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gZ_ >> gWL_ >> _couplingValueLeftIm >> _couplingValueLeftRe >> _couplingValueRightIm >> _couplingValueRightRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<HalfHalfOneEWSplitFn,SplittingFunction>
describeHerwigHalfHalfOneEWSplitFn("Herwig::HalfHalfOneEWSplitFn", "HwShower.so");

void HalfHalfOneEWSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneEWSplitFn> documentation
    ("The HalfHalfOneEWSplitFn class implements the splitting q->qW and q->qZ");

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueLeftIm
    ("CouplingValue.Left.Im",
     "The numerical value (imaginary part) of the left-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueLeftIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueLeftRe
    ("CouplingValue.Left.Re",
     "The numerical value (real part) of the left-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueLeftRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueRightIm
    ("CouplingValue.Right.Im",
     "The numerical value (imaginary part) of the right-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueRightIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<HalfHalfOneEWSplitFn, double> interfaceCouplingValueRightRe
    ("CouplingValue.Right.Re",
     "The numerical value (real part) of the right-handed splitting coupling to be imported for BSM splittings",
     &HalfHalfOneEWSplitFn::_couplingValueRightRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}

void HalfHalfOneEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // left-handled W coupling
  gWL_ = 1./sqrt(2.*sw2);
  // Z couplings
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    gZ_[2*ix-1]  = make_pair(fact*(sm->vd()  + sm->ad()),
			     fact*(sm->vd()  - sm->ad() ));
    gZ_[2*ix  ]  = make_pair(fact*(sm->vu()  + sm->au() ),
			     fact*(sm->vu()  - sm->au() ));
    gZ_[2*ix+9 ] = make_pair(fact*(sm->ve()  + sm->ae() ),
			     fact*(sm->ve()  - sm->ae() ));
    gZ_[2*ix+10] = make_pair(fact*(sm->vnu() + sm->anu()),
			     fact*(sm->vnu() - sm->anu()));
  }
}

void HalfHalfOneEWSplitFn::getCouplings(Complex & gL, Complex & gR, const IdList & ids) const {
  if(_couplingValueLeftRe!=0||_couplingValueLeftIm!=0||_couplingValueRightRe!=0||_couplingValueRightIm!=0) {
    double e = sqrt(4.*Constants::pi*generator()->standardModel()
              ->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
    // a factor e is factored out since its already accounted for
    gL = Complex(_couplingValueLeftRe,_couplingValueLeftIm)/e;
    gR = Complex(_couplingValueRightRe,_couplingValueRightIm)/e;
  }
  else if(ids[2]->id()==ParticleID::Z0) {
    map<long,pair<double,double> >::const_iterator it = gZ_.find(abs(ids[0]->id()));
    assert(it!=gZ_.end());
    gL = it->second.first ;
    gR = it->second.second;
  }
  else if(abs(ids[2]->id())==ParticleID::Wplus) {
    gL = gWL_;
  }
  else
    assert(false);
  if(ids[0]->id()<0) swap(gL,gR);
}

double HalfHalfOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  double val = (norm(gL)*abs(rho(0,0))+norm(gR)*abs(rho(1,1)))*(1.+sqr(z))/(1.-z);
  Energy m0, m1, m2;
  if(mass) {
    m0 = ids[0]->mass();
    m1 = ids[1]->mass();
    m2 = ids[2]->mass();
    double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);
    val += (norm(gL)*abs(rho(0,0))+norm(gR)*abs(rho(1,1)))*((1.+sqr(z))/(1.-z)*sqr(m0t)-(1.+z)/(1.-z)*sqr(m1t)-sqr(m2t))
           + (norm(gR)*abs(rho(0,0))+norm(gL)*abs(rho(1,1)))*z*sqr(m0t)
           - 2.*(gR*conj(gL)).real()*(abs(rho(1,1))+abs(rho(0,0)))*m0t*m1t;
  }
  return colourFactor(ids)*val;
}


double HalfHalfOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  return 2.*max(norm(gL),norm(gR))*colourFactor(ids)/(1.-z); //FIXME//
}

double HalfHalfOneEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  double val = (norm(gL)*abs(rho(0,0))+norm(gR)*abs(rho(1,1)))*(1.+sqr(z));
  Energy m0, m1, m2;
  if(mass) {
    m0 = ids[0]->mass();
    m1 = ids[1]->mass();
    m2 = ids[2]->mass();
    double m0t = m0/sqrt(t), m1t = m1/sqrt(t), m2t = m2/sqrt(t);
    val += (norm(gL)*abs(rho(0,0))+norm(gR)*abs(rho(1,1)))*((1+sqr(z))*sqr(m0t)-(1.+z)*sqr(m1t)-(1.-z)*sqr(m2t))
           + (norm(gR)*abs(rho(0,0))+norm(gL)*abs(rho(1,1)))*z*(1.-z)*sqr(m0t)
           - 2.*(gR*conj(gL)).real()*(abs(rho(1,1))+abs(rho(0,0)))*(1.-z)*m0t*m1t;
  }
  val /= max(norm(gR),norm(gL));
  return 0.5*val;
}

double HalfHalfOneEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  double pre = colourFactor(ids)*max(norm(gL),norm(gR));
  switch (PDFfactor) {
  case 0:
    return -2.*pre*Math::log1m(z);
  case 1:
    return  2.*pre*log(z/(1.-z));
  case 2:
    return  2.*pre/(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfHalfOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  double pre = colourFactor(ids)*max(norm(gL),norm(gR));
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/pre);
  case 1:
    return 1./(1.-exp(-0.5*r/pre));
  case 2:
    return 1.-2.*pre/r;
  case 3:
  default:
    throw Exception() << "HalfHalfOneEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfHalfOneEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[2]->iSpin()==PDT::Spin1 && !(_couplingValueLeftRe==0 && _couplingValueLeftIm==0 && _couplingValueRightRe==0 && _couplingValueRightIm==0)) {
    if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge()) return false;
    if((abs(ids[0]->id())>=1 && abs(ids[0]->id())<=6) && (abs(ids[1]->id())>=1 && abs(ids[1]->id())<=6)) return true;
  }
  else if(ids[2]->id()==ParticleID::Z0) {
    if(ids[0]->id()==ids[1]->id() &&
       ((ids[0]->id()>=1 && ids[0]->id()<=6) || (ids[0]->id()>=11&&ids[0]->id()<=16) )) return true;
  }
  else if(abs(ids[2]->id())==ParticleID::Wplus) {
    if(!((ids[0]->id()>=1 && ids[0]->id()<=6) || (ids[0]->id()>=11&&ids[0]->id()<=16) )) return false;
    if(!((ids[1]->id()>=1 && ids[1]->id()<=6) || (ids[1]->id()>=11&&ids[1]->id()<=16) )) return false;
    if(ids[0]->id()+1!=ids[1]->id() && ids[0]->id()-1!=ids[1]->id()) return false;
    int out = ids[1]->iCharge()+ids[2]->iCharge();
    if(ids[0]->iCharge()==out) return true;
  }
  return false;
}

vector<pair<int, Complex> >
HalfHalfOneEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

vector<pair<int, Complex> >
HalfHalfOneEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

DecayMEPtr HalfHalfOneEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Energy m0 = ids[0]->mass();
  Energy m1 = ids[1]->mass();
  Energy m2 = ids[2]->mass();
  Complex gL(0.,0.),gR(0.,0.);
  getCouplings(gL,gR,ids);
  Energy den = sqrt(2.)*sqrt(t);
  Energy pt  = sqrt(z*(1.-z)*t+z*(1.-z)*sqr(m0)-(1.-z)*sqr(m1)-z*sqr(m2));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  (*kernal)(1,1,2) = sqrt(2.)*gR*pt/sqrt(z)/(1.-z)*cphase/den;
  (*kernal)(1,1,1) = -2.*gR*sqrt(z)*m2/(1.-z)/den;
  (*kernal)(1,1,0) = -sqrt(2.*z)*gR*pt/(1.-z)*phase/den;
  (*kernal)(1,0,2) = -sqrt(2.)*(gL*z*m0-gR*m1)/sqrt(z)/den;
  (*kernal)(1,0,1) = 0;
  (*kernal)(1,0,0) = 0;
  (*kernal)(0,1,2) = 0;
  (*kernal)(0,1,1) = 0;
  (*kernal)(0,1,0) = -sqrt(2.)*(gR*z*m0-gL*m1)/sqrt(z)/den;
  (*kernal)(0,0,2) = sqrt(2.*z)*gL*pt/(1.-z)*cphase/den;
  (*kernal)(0,0,1) = -2.*gL*sqrt(z)*m2/(1.-z)/den;
  (*kernal)(0,0,0) = -sqrt(2.)*gL*pt/sqrt(z)/(1.-z)*phase/den;

  // return the answer
  return kernal;
}
