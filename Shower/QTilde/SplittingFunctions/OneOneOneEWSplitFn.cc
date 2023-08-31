// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneOneEWSplitFn class.
//

#include "OneOneOneEWSplitFn.h"
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

IBPtr OneOneOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneOneOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneOneOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gWWG_ << gWWZ_ << _theSM << _couplingValueIm << _couplingValueRe;
}

void OneOneOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gWWG_ >> gWWZ_ >> _theSM >> _couplingValueIm >> _couplingValueRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneOneOneEWSplitFn,SplittingFunction>
describeHerwigOneOneOneEWSplitFn("Herwig::OneOneOneEWSplitFn", "HwShower.so");


void OneOneOneEWSplitFn::Init() {

  static ClassDocumentation<OneOneOneEWSplitFn> documentation
    ("The OneOneOneEWSplitFn class implements the splitting W->WG, W->WZ and Z->ZZ");

  static Parameter<OneOneOneEWSplitFn, double> interfaceCouplingValueIm
    ("CouplingValue.Im",
     "The numerical value (imaginary part) of the splitting coupling to be imported for BSM splittings",
     &OneOneOneEWSplitFn::_couplingValueIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<OneOneOneEWSplitFn, double> interfaceCouplingValueRe
    ("CouplingValue.Re",
     "The numerical value (real part) of the splitting coupling to be imported for BSM splittings",
     &OneOneOneEWSplitFn::_couplingValueRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}


void OneOneOneEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // WWZ coupling
  gWWZ_ = sqrt((1.-sw2)/sw2);
  // WWG coupling
  gWWG_ = 1.;
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}


void OneOneOneEWSplitFn::getCouplings(Complex & gvvv, const IdList & ids) const {
  if(_couplingValueIm!=0||_couplingValueRe!=0) {
    double e  = sqrt(4.*Constants::pi*generator()->standardModel()
              ->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
    gvvv = Complex(_couplingValueRe,_couplingValueIm)/e;
  }
  // Z > WW
  else if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                       && abs(ids[2]->id())==ParticleID::Wplus){
    gvvv = Complex(0.,gWWZ_);
  }
  // W > WG
  else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                               && ids[2]->id()==ParticleID::gamma){
    gvvv = Complex(0.,gWWG_);
  }
  // W > WZ
  else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                               && ids[2]->id()==ParticleID::Z0){
    gvvv = Complex(0.,gWWZ_);
  }
  else
    assert(false);
}


double OneOneOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  Complex gvvv(0.,0.);
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
    val += 4.*m0t2*sqr(1.-z)*abs_rho_11 + (2*(m0t2*sqr(1.-(1.-z)*z)
        -m2t2*(1.-sqr(1.-z)*z) -m1t2*(1.-(1.-z)*sqr(z)))*(abs_rho_00+abs_rho_22))
        /((1.-z)*z);
  }
  return norm(gvvv)*val;
}


double OneOneOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex gvvv(0.,0.);
  getCouplings(gvvv,ids);
  double val = norm(gvvv)*(2./(z*(1.-z)));
  return val;
}


double OneOneOneEWSplitFn::ratioP(const double z, const Energy2 t,
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
    val += (4.*m0t2*sqr(1.-z)*abs_rho_11 + (2*(m0t2*sqr(1.-(1.-z)*z)
        -m2t2*(1.-sqr(1.-z)*z) -m1t2*(1.-(1.-z)*sqr(z)))*(abs_rho_00+abs_rho_22))
        /((1.-z)*z))/(2./(z*(1.-z)));
  }
  return val;
}


double OneOneOneEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex gvvv(0.);
  getCouplings(gvvv,ids);
  double pre = norm(gvvv);
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
    throw Exception() << "OneOneOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneOneOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex gvvv(0.);
  getCouplings(gvvv,ids);
  double pre = norm(gvvv);
  switch (PDFfactor) {
  case 0:
    return exp(0.5*r/pre)/(1.+exp(0.5*r/pre));
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneOneEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  // Z > WW
  if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                       && ids[1]->id()==-ids[2]->id())
    return true;

  if(abs(ids[0]->id())==ParticleID::Wplus) {
    // W > WG
    if(ids[1]->id()==ids[0]->id() && ids[2]->id()==ParticleID::gamma)
      return true;
    // W > WZ
    if(ids[1]->id()==ids[0]->id() && ids[2]->id()==ParticleID::Z0)
      return true;
  }
  return false;
}


vector<pair<int, Complex> >
OneOneOneEWSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


vector<pair<int, Complex> >
OneOneOneEWSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}


DecayMEPtr OneOneOneEWSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Complex gvvv(0.);
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
