// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneOneQEDSplitFn class.
//

#include "OneOneOneQEDSplitFn.h"
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

IBPtr OneOneOneQEDSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneOneOneQEDSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneOneOneQEDSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gWWG_ << _theSM;
}

void OneOneOneQEDSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gWWG_ >> _theSM;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneOneOneQEDSplitFn,Sudakov1to2FormFactor>
describeHerwigOneOneOneQEDSplitFn("Herwig::OneOneOneQEDSplitFn", "HwShower.so");


void OneOneOneQEDSplitFn::Init() {

  static ClassDocumentation<OneOneOneQEDSplitFn> documentation
    ("The OneOneOneQEDSplitFn class implements the gamma->WW EW splitting.");

}


void OneOneOneQEDSplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // WWG coupling
  gWWG_ = 1.;
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}


void OneOneOneQEDSplitFn::getCouplings(double & gvvv, const IdList & ids) const {
  // G > WW
  if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                     && abs(ids[2]->id())==ParticleID::Wplus){
    gvvv = gWWG_;
  }
  else
    assert(false);
}


double OneOneOneQEDSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  double abs_rho_00 = abs(rho(0,0));
  double abs_rho_11 = abs(rho(1,1));
  double abs_rho_22 = abs(rho(2,2));
  // massless limit
  double val = ((2.*sqr(1.-(1.-z)*z))/((1.-z)*z))*(abs_rho_00+abs_rho_22);
  // massive limits
  if(mass) {
    double m0t2 = sqr(ids[0]->mass())/t;
    double m1t2 = sqr(ids[1]->mass())/t;
    double m2t2 = sqr(ids[2]->mass())/t;
    val += (-2.*(m2t2*(1.-sqr(1.-z)*z)+m1t2*(1.-(1.-z)*sqr(z)))*(abs_rho_00+abs_rho_22))/((1.-z)*z)
      + (2.*m0t2*(2.*pow(1.-z,3)*z*abs_rho_11+sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22)))/((1.-z)*z);
  }
  return sqr(gvvv)*val;
}


double OneOneOneQEDSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  return sqr(gvvv)*(2./(z*(1.-z)));
}


double OneOneOneQEDSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double val(0.);
  double abs_rho_00 = abs(rho(0,0));
  double abs_rho_11 = abs(rho(1,1));
  double abs_rho_22 = abs(rho(2,2));
  // massless limit
  val = sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22);
  // massive limit
  if(mass) {
        double m0t2 = sqr(ids[0]->mass())/t;
        double m1t2 = sqr(ids[1]->mass())/t;
        double m2t2 = sqr(ids[2]->mass())/t;
    val += -(m2t2*(1.-sqr(1.-z)*z) + m1t2*(1.-(1.-z)*sqr(z)))*(abs_rho_00+abs_rho_22)
         + m0t2*(2.*pow(1.-z,3)*z*abs_rho_11+sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22));
  }
  return val;
}


double OneOneOneQEDSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  double pre = sqr(gvvv);
  switch (PDFfactor) {
  case 0:
    return 2.*pre*(log(z)-log(1.-z));
  case 1:
    //return -2.*pre*(1./z+log(1.-z)-log(z));
  case 2:
    //return 2.*pre*(2.*log(z)+(2.*z-1.)/(z*(1.-z))-2.*log(1.-z));
  case 3:
    //return 2.*pre*(1./(1.-z)-1./z-2.*log(1.-z)+2.*log(z));
  default:
    throw Exception() << "OneOneOneQEDSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneOneOneQEDSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  double pre = sqr(gvvv);
  switch (PDFfactor) {
  case 0:
    return exp(0.5*r/pre)/(1.+exp(0.5*r/pre));
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneQEDSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneOneQEDSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                     && ids[1]->id()==-ids[2]->id())
    return true;
  return false;
}


vector<pair<int, Complex> >
OneOneOneQEDSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


vector<pair<int, Complex> >
OneOneOneQEDSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


DecayMEPtr OneOneOneQEDSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  double gvvv(0.);
  getCouplings(gvvv,ids);
  // defining dummies
  double m0t = ids[0]->mass()/sqrt(t);
  double m1t = ids[1]->mass()/sqrt(t);
  double m2t = ids[2]->mass()/sqrt(t);
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);

  double z1_z = z*(1.-z);
  double sqrtmass = sqrt(sqr(m0t)-sqr(m1t)/z-sqr(m2t)/(1.-z)+1.);
  double r2   = sqrt(2.);
  // assign kernel
  (*kernal)(0,0,0) = gvvv*phase*(1./sqrt(z1_z))*sqrtmass;
  (*kernal)(0,0,1) = gvvv*r2*m2t*(z/(1.-z)); //2>4
  (*kernal)(0,0,2) = -gvvv*cphase*sqrt(z/(1.-z))*sqrtmass;
  (*kernal)(0,1,0) = -gvvv*r2*m1t*(1.-z)/z; //2>4
  (*kernal)(0,1,1) = 0.;
  (*kernal)(0,1,2) = 0.;
  (*kernal)(0,2,0) = -gvvv*(1.-z)*cphase*sqrt((1.-z)/z)*sqrtmass;
  (*kernal)(0,2,1) = 0.;
  (*kernal)(0,2,2) = 0.;

  (*kernal)(1,0,0) = 0.;
  (*kernal)(1,0,1) = 0.; //2>4
  (*kernal)(1,0,2) = -gvvv*r2*m0t*(1.-z); //2>4
  (*kernal)(1,1,0) = 0.; //221>441
  (*kernal)(1,1,1) = 0.; //222>444
  (*kernal)(1,1,2) = 0.; //223>443
  (*kernal)(1,2,0) = -gvvv*r2*m0t*(1.-z); //2>4
  (*kernal)(1,2,1) = 0.; //2>4
  (*kernal)(1,2,2) = 0.; //2>4

  (*kernal)(2,0,0) = 0.;
  (*kernal)(2,0,1) = 0.;
  (*kernal)(2,0,2) = gvvv*(1.-z)*phase*sqrt((1.-z)/z)*sqrtmass;
  (*kernal)(2,1,0) = 0.;
  (*kernal)(2,1,1) = 0.; //2>4
  (*kernal)(2,1,2) = -gvvv*r2*m1t*((1.-z)/z);//2>4
  (*kernal)(2,2,0) = gvvv*phase*sqrt(z/(1.-z))*sqrtmass;
  (*kernal)(2,2,1) = gvvv*r2*m2t*(z/(1.-z)); //2>4
  (*kernal)(2,2,2) = -gvvv*cphase*(1./sqrt(z1_z))*sqrtmass;

  // return the answer
  return kernal;
}
