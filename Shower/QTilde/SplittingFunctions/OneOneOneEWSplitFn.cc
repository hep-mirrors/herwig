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

using namespace Herwig;

IBPtr OneOneOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneOneOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneOneOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gWWG_ << gWWZ_ << _theSM;
}

void OneOneOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gWWG_ >> gWWZ_ >> _theSM;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneOneOneEWSplitFn,SplittingFunction>
describeHerwigOneOneOneEWSplitFn("Herwig::OneOneOneEWSplitFn", "HwShower.so");


void OneOneOneEWSplitFn::Init() {

  static ClassDocumentation<OneOneOneEWSplitFn> documentation
    ("The OneOneOneEWSplitFn class implements the splitting W->WG, W->WZ and Z->ZZ");

}


void OneOneOneEWSplitFn::doinit() {
  SplittingFunction::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // WWZ coupling
  gWWZ_ = 1.;
  // WWG coupling
  gWWG_ = sqrt((1.-sw2)/sw2);
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}


void OneOneOneEWSplitFn::getCouplings(double & gvvv, const IdList & ids) const {
  // G > WW
  if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                     && abs(ids[2]->id())==ParticleID::Wplus){
    gvvv = gWWG_;
  }
  // Z > WW
  else if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                       && abs(ids[2]->id())==ParticleID::Wplus){
    gvvv = gWWZ_;
  }
  // W > WG
  else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                               && ids[2]->id()==ParticleID::gamma){
    gvvv = gWWG_;
  }
  // W > WZ
  else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                               && ids[2]->id()==ParticleID::Z0){
    gvvv = gWWZ_;
  }
  else
    assert(false);
}


double OneOneOneEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  double val(0.);
  double abs_rho_00 = sqrt(norm(rho(0,0)));
  double abs_rho_11 = sqrt(norm(rho(1,1)));
  double abs_rho_22 = sqrt(norm(rho(2,2)));
  val = ((2.*sqr(1.-(1.-z)*z))/((1.-z)*z))*(abs_rho_00+abs_rho_22);
  if(mass) {
    double mWt2 = sqr(getParticleData(ParticleID::Wplus)->mass())/t;
    double mZt2 = sqr(getParticleData(ParticleID::Z0)->mass())/t;
    // G > WW
    if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                       && abs(ids[2]->id())==ParticleID::Wplus){
      val += (-2.*mWt2*(2.-(1.-z)*z)*(abs_rho_00+abs_rho_22))/((1.-z)*z);
    }
    // Z > WW
    else if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                         && abs(ids[2]->id())==ParticleID::Wplus){
      val += (-2.*mWt2*(2.-(1.-z)*z)*(abs_rho_00+abs_rho_22))/((1.-z)*z);
      val += (-2.*mZt2*(z*(-2.+z*(8.+z*(-13.-2.*(-4.+z)*z)))*abs_rho_11
           - (1.-z)*sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22)))/(sqr(1.-z)*z);
    }
    // W > WG
    else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                                 && ids[2]->id()==ParticleID::gamma){
      val += 4.*mWt2*sqr(1.-z)*abs_rho_11;
      val -= 2.*mWt2*(1.+sqr(1.-z))*(abs_rho_00 + abs_rho_22);
    }
    // W > WZ
    else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                                 && ids[2]->id()==ParticleID::Z0){
      val += (2.*mZt2*(pow(z,3)*abs_rho_11 + (1.-z)*(-1.+sqr(1.-z)*z)
           * (abs_rho_00+abs_rho_22)))/(sqr(1.-z)*z);
      val += -2.*mWt2*(-2.*sqr(1.-z)*abs_rho_11+(2.+(-2.+z)*z)
           * (abs_rho_00+abs_rho_22));
    }
  }
  return sqr(gvvv)*val;
}


double OneOneOneEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  return sqr(gvvv)*(2.*sqr(gvvv)/(z*(1.-z)));
}


double OneOneOneEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double gvvv(0.);
  getCouplings(gvvv,ids);
  double abs_rho_00 = sqrt(norm(rho(0,0)));
  double abs_rho_11 = sqrt(norm(rho(1,1)));
  double abs_rho_22 = sqrt(norm(rho(2,2)));
  double val = sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22);
  if(mass) {
    double mWt2 = sqr(getParticleData(ParticleID::Wplus)->mass())/t;
    double mZt2 = sqr(getParticleData(ParticleID::Z0)->mass())/t;
    // G > WW
    if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                       && abs(ids[2]->id())==ParticleID::Wplus) {
      val += (-2.*mWt2*(2.-(1.-z)*z)*(abs_rho_00+abs_rho_22))/((1.-z)*z);
    }
    // Z > WW
    else if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                         && abs(ids[2]->id())==ParticleID::Wplus){
      val += (-2.*mWt2*(2.-(1.-z)*z)*(abs_rho_00+abs_rho_22))/((1.-z)*z);
      val += (-2.*mZt2*(z*(-2.+z*(8.+z*(-13.-2.*(-4.+z)*z)))*abs_rho_11
           - (1.-z)*sqr(1.-(1.-z)*z)*(abs_rho_00+abs_rho_22)))/(sqr(1.-z)*z);
    }
    // W > WG
    else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                                 && ids[2]->id()==ParticleID::gamma) {
      val += 4.*mWt2*sqr(1.-z)*abs_rho_11;
      val -= 2.*mWt2*(1.+sqr(1.-z))*(abs_rho_00 + abs_rho_22);
    }
    // W > WZ
    else if(abs(ids[0]->id())==ParticleID::Wplus && abs(ids[1]->id())==ParticleID::Wplus
                                                 && ids[2]->id()==ParticleID::Z0) {
      val += (2.*mZt2*(pow(z,3)*abs_rho_11 + (1.-z)*(-1.+sqr(1.-z)*z)
             * (abs_rho_00+abs_rho_22)))/(sqr(1.-z)*z);
      val += -2.*mWt2*(-2.*sqr(1.-z)*abs_rho_11+(2.+(-2.+z)*z)
             * (abs_rho_00+abs_rho_22));
    }
  }
  return sqr(gvvv)*val;
}


double OneOneOneEWSplitFn::integOverP(const double z,
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
    throw Exception() << "OneOneOneEWSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneOneOneEWSplitFn::invIntegOverP(const double r, const IdList & ids,
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
    throw Exception() << "OneOneOneEWSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneOneEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]->id()==ParticleID::gamma && abs(ids[1]->id())==ParticleID::Wplus
                                     && ids[1]->id()==-ids[2]->id())
    return true;

  if(ids[0]->id()==ParticleID::Z0 && abs(ids[1]->id())==ParticleID::Wplus
                                       && ids[1]->id()==-ids[2]->id())
    return true;

  if(abs(ids[0]->id())==ParticleID::Wplus) {
    if(ids[1]->id()==ids[0]->id() && ids[2]->id()==ParticleID::gamma)
      return true;
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
  (*kernal)(0,0,0) = gvvv*(phase/sqrt(z1_z))*sqrtmass;
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
  (*kernal)(1,1,0) = gvvv*phase*(m0t/m1t)*sqrt(z/(1.-z));   //221>421
  (*kernal)(1,1,1) = gvvv*r2*(m0t*m2t/m1t)*(z/(1.-z));      //222>424
  (*kernal)(1,1,2) = -gvvv*cphase*(m0t/m1t)*sqrt(z/(1.-z)); //223>423
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
  (*kernal)(2,2,2) = -gvvv*(cphase/sqrt(z1_z))*sqrtmass;

  // return the answer
  return kernal;
}
