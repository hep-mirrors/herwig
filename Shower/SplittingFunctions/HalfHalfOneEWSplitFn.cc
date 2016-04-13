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

using namespace Herwig;

IBPtr HalfHalfOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HalfHalfOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void HalfHalfOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gZ_ << gWL_;
}

void HalfHalfOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gZ_ >> gWL_;
}

// The following static variable is needed for the type description system in ThePEG. 
DescribeClass<HalfHalfOneEWSplitFn,SplittingFunction>
describeHerwigHalfHalfOneEWSplitFn("Herwig::HalfHalfOneEWSplitFn", "HwShower.so");

void HalfHalfOneEWSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneEWSplitFn> documentation
    ("The HalfHalfOneEWSplitFn class implements the splitting q->qWand q->qZ");

}

void HalfHalfOneEWSplitFn::doinit() {
  cerr << "testing in do init\n";
  SplittingFunction::doinit();
  cerr << "testing in do init\n";
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

void HalfHalfOneEWSplitFn::getCouplings(double & gL, double & gR, const IdList & ids) const {
  if(ids[2]==ParticleID::Z0) {
    map<long,pair<double,double> >::const_iterator it = gZ_.find(abs(ids[0]));
    assert(it!=gZ_.end());
    gL = it->second.first ;
    gR = it->second.second;
  }
  else if(abs(ids[2])==ParticleID::Wplus) {
    gL = gWL_;
  }
  else
    assert(false);
}

double HalfHalfOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass) const {
  double gL(0.),gR(0.);
  getCouplings(gL,gR,ids);
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[2])->mass();  
    val -= sqr(m)/t;
  }
  val *= (sqr(gL)+sqr(gR));
  return colourFactor(ids)*val;
}


double HalfHalfOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gL(0.),gR(0.);
  getCouplings(gL,gR,ids);
  return 2.*sqr(max(gL,gR))*colourFactor(ids)/(1.-z); 
}

double HalfHalfOneEWSplitFn::ratioP(const double z, const Energy2 t,
				  const IdList & ids, const bool mass) const {
  double gL(0.),gR(0.);
  getCouplings(gL,gR,ids);
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[2])->mass();  
    val -= (1.-z)*sqr(m)/t;
  }
  val *= 0.5*(sqr(gL)+sqr(gR))/sqr(max(gL,gR));
  return val;
} 

double HalfHalfOneEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gL(0.),gR(0.);
  getCouplings(gL,gR,ids);
  double pre = colourFactor(ids)*sqr(max(gL,gR));
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
  double gL(0.),gR(0.);
  getCouplings(gL,gR,ids);
  double pre = colourFactor(ids)*sqr(max(gL,gR));
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
  if(ids[2]==ParticleID::Z0) {
    if(ids[0]==ids[1] && 
       ((ids[0]>=1 && ids[0]<=6) || (ids[0]>=11&&ids[0]<=16) )) return true;
  }
  else if(abs(ids[2])==ParticleID::Wplus) {
    if(!((ids[0]>=1 && ids[0]<=6) || (ids[0]>=11&&ids[0]<=16) )) return false;
    if(!((ids[1]>=1 && ids[1]<=6) || (ids[1]>=11&&ids[1]<=16) )) return false;
    if(ids[0]+1!=ids[1] && ids[0]-1!=ids[1]) return false;
    int out = getParticleData(ids[1])->iCharge()+getParticleData(ids[2])->iCharge();
    if(getParticleData(ids[0])->iCharge()==out) return true;
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
                                             bool timeLike) {
//   // calculate the kernal
//   DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
//   Energy m = !timeLike ? ZERO : getParticleData(ids[0])->mass();
//   double mt = m/sqrt(t);
//   double root = sqrt(1.-(1.-z)*sqr(m)/z/t);
//   double romz = sqrt(1.-z); 
//   double rz   = sqrt(z);
//   Complex phase = exp(Complex(0.,1.)*phi);
//   (*kernal)(0,0,0) = -root/romz*phase;
//   (*kernal)(1,1,2) =  -conj((*kernal)(0,0,0));
//   (*kernal)(0,0,2) =  root/romz*z/phase;
//   (*kernal)(1,1,0) = -conj((*kernal)(0,0,2));
//   (*kernal)(1,0,2) =  mt*(1.-z)/rz;
//   (*kernal)(0,1,0) =  conj((*kernal)(1,0,2));
//   (*kernal)(0,1,2) =  0.;
//   (*kernal)(1,0,0) =  0.;
//   return kernal;
}
