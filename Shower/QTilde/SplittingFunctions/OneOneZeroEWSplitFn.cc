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
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr OneOneZeroEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr OneOneZeroEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void OneOneZeroEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << gWWH_ << gZZH_ << _theSM << _couplingValueIm << _couplingValueRe;
}

void OneOneZeroEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> gWWH_ >> gZZH_ >> _theSM >> _couplingValueIm >> _couplingValueRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<OneOneZeroEWSplitFn,Sudakov1to2FormFactor>
describeHerwigOneOneZeroEWSplitFn("Herwig::OneOneZeroEWSplitFn", "HwShower.so");


void OneOneZeroEWSplitFn::Init() {

  static ClassDocumentation<OneOneZeroEWSplitFn> documentation
    ("The OneOneZeroEWSplitFn class implements the splittings W->WH and Z->ZH");

  static Parameter<OneOneZeroEWSplitFn, double> interfaceCouplingValueIm
    ("CouplingValue.Im",
     "The numerical value (imaginary part) of the splitting coupling to be imported for BSM splittings",
     &OneOneZeroEWSplitFn::_couplingValueIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<OneOneZeroEWSplitFn, double> interfaceCouplingValueRe
    ("CouplingValue.Re",
     "The numerical value (real part) of the splitting coupling to be imported for BSM splittings",
     &OneOneZeroEWSplitFn::_couplingValueRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}


void OneOneZeroEWSplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sm->sin2ThetaW();
  // W -> W H coupling      g = e/sin theta_W
  gWWH_ = 1./sqrt(sw2);
  // Z -> Z H coupling
  gZZH_ = 1./(sqrt(sw2)*sqrt(1.-sw2));
  // to employ running masses, wherever needed
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
}

void OneOneZeroEWSplitFn::getCouplings(Complex & g, const IdList & ids) const {
  if(_couplingValueIm!=0||_couplingValueRe!=0) {
    sm_ = false;
    double e  = sqrt(4.*Constants::pi*generator()->standardModel()
              ->alphaEM(sqr(getParticleData(ParticleID::Z0)->mass())));
    g = Complex(_couplingValueRe,_couplingValueIm)/e;
  }
  else {
    if(abs(ids[0]->id())==ParticleID::Wplus)
      g = Complex(0.,gWWH_);
    else if(ids[0]->id()==ParticleID::Z0)
      g = Complex(0.,gZZH_);
    else
      assert(false);
  }
}

double OneOneZeroEWSplitFn::P(const double z, const Energy2 t,
			       const IdList &ids, const bool mass, const RhoDMatrix & rho) const {
  Complex gvvh(0.,0.);
  getCouplings(gvvh,ids);
  double rho00 = abs(rho(0,0)), rho11 = abs(rho(1,1)), rho22 = abs(rho(2,2));
  double val;
  if(sm_) { 
    val = (1.-z)*(2.*rho11+sqr(z)*(rho00+rho22));
    // the massive limit
    if( mass ) {
      // get the running mass
      double mVt = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
      double mHt = _theSM->mass(t,getParticleData(ids[2]->id()))/sqrt(t);
      val -= sqr(mHt)*(2.*rho11+sqr(z)*(rho00+rho22))
           + sqr(mVt)/z*(rho11*(2.*sqr(z)-4.*z+2.)+(rho00+rho22)*sqr(z)*(sqr(z)-2.*z-1.));
    }
    val *= norm(gvvh)/(4.*z);
  }
  else {
    double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
    double m0t = ids[0]->mass()/sqrt(t), m1t = ids[1]->mass()/sqrt(t), m2t = ids[2]->mass()/sqrt(t);
    if(mass) {
      if(m0t!=0.&&m1t!=0.)
        val = (norm(gvvh)/(4.*sqr(m1)))*(rho00+rho22)*(z*(1.-z)*(1.+sqr(m0t))+(1.+z)*sqr(m1t)-z*sqr(m2t))
            + (norm(gvvh)/(2.*sqr(m0)))*rho11/sqr(z)*(z*(1.-z)*(1.+sqr(m0t))-(1.-z)*sqr(m1t)-z*sqr(m2t));
      else if(m1t!=0.)
        val = (norm(gvvh)/(4.*sqr(m1)))*(z*(1.-z)+(1.+z)*sqr(m1t)-z*sqr(m2t));
      else if(m0t!=0.)
        val = (norm(gvvh)/(2.*sqr(m0)))*((rho00+rho22)*sqr(m0t)
            + rho11/sqr(z)*(z*(1-z)+z*(1.-z)*sqr(m0t)-z*sqr(m2t)));
      else
        val = norm(gvvh)/(2.*t)*GeV2;
    }
    else {
      val = norm(gvvh)/(2.*t)*GeV2;
    }
  }
  return val;
}


double OneOneZeroEWSplitFn::overestimateP(const double z,
					   const IdList & ids) const {
  Complex gvvh(0.,0.);
  getCouplings(gvvh,ids);
  if(sm_)
    return norm(gvvh)/(2.*z);
  else {
    double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
    if(m0!=0.&&m1!=0.) {
      if(m0<m1) gvvh /= m0;
      else gvvh /= m1;
    }
    else if(m1!=0.) gvvh /= m1;
    else if(m0!=0.) gvvh /= m0;
    return norm(gvvh)/(2.*z*(1.-z));
  }
}


double OneOneZeroEWSplitFn::ratioP(const double z, const Energy2 t,
				    const IdList & ids, const bool mass,
				    const RhoDMatrix & rho) const {
  double rho00 = abs(rho(0,0));
  double rho11 = abs(rho(1,1));
  double rho22 = abs(rho(2,2));
  double val;
  if(sm_) {
    // ratio in the massless limit
    val = (1.-z)*(2.*rho11+sqr(z)*(rho00+rho22))/2.;
    // the massive limit
    if(mass){
      // get the running mass
      double mVt = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
      double mHt = _theSM->mass(t,getParticleData(ids[2]->id()))/sqrt(t);
      val -= sqr(mHt)*(2.*rho11+sqr(z)*(rho00+rho22))
           + sqr(mVt)/z*(rho11*(2.*sqr(z)-4.*z+2.)+(rho00+rho22)*sqr(z)*(sqr(z)-2.*z-1.));
    }
    val /= 2.;
  }
  else {
    if(mass) {
      double m0t = ids[0]->mass()/sqrt(t), m1t = ids[1]->mass()/sqrt(t), m2t = ids[2]->mass()/sqrt(t);
      if(m0t!=0.&&m1t!=0.) {
        double mR1, mR2; // square of mass ratio between m0 and m1
        if(m0t<m1t) {
          mR1 = sqr(m0t)/sqr(m1t); mR2 = 1.;
        }
        else {
          mR1 = 1.; mR2 = sqr(m1t)/sqr(m0t);
        }
        val = mR1*(rho00+rho22)/2.*(z*(1.-z)*(1.+sqr(m0t))+(1.+z)*sqr(m1t)-z*sqr(m2t))
            + mR2*rho11/sqr(z)*(z*(1.-z)*(1.+sqr(m0t))-(1.-z)*sqr(m1t)-z*sqr(m2t));
      }
      else if(m1t!=0.)
        val = (z*(1.-z)+(1.+z)*sqr(m1t)-z*sqr(m2t))/2.;
      else if(m0t!=0.)
        val = (rho00+rho22)*sqr(m0t)+2.*rho11*(z*(1-z)+z*(1.-z)*sqr(m0t)-z*sqr(m2t))/sqr(z);
      else
        val = 1./t*GeV2;
    }
    else {
      double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
      val = 1./t*GeV2;
      if(m0!=0.&&m1!=0.) {
        if(m0<m1) val *= sqr(m0);
        else val *= sqr(m1);
      }
      else if(m1!=0.) val *= sqr(m1);
      else if(m0!=0.) val *= sqr(m0);
    }
    val *= z*(1.-z);
  }
  return val;
}


double OneOneZeroEWSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  Complex gvvh(0.,0.);
  getCouplings(gvvh,ids);
  double pre = norm(gvvh);
  if(sm_) {
    switch (PDFfactor) {
    case 0:
      return pre*log(z)/2.;
    case 1:
    case 2:
    case 3:
    default:
      throw Exception() << "OneOneZeroEWSplitFn::integOverP() invalid PDFfactor = "
  		      << PDFfactor << Exception::runerror;
    }
  }
  else {
    double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
    if(m0!=0.&&m1!=0.) {
      if(m0<m1) pre /= sqr(m0);
      else pre /= sqr(m1);
    }
    else if(m1!=0.) pre /= sqr(m1);
    else if(m0!=0.) pre /= sqr(m0);
    switch (PDFfactor) {
    case 0:
      return pre*(log(z)-log(1.-z))/2.;
    case 1:
    case 2:
    case 3:
    default:
      throw Exception() << "OneOneZeroEWSplitFn::integOverP() invalid PDFfactor = "
  		      << PDFfactor << Exception::runerror;
    }
  }
}

double OneOneZeroEWSplitFn::invIntegOverP(const double r, const IdList & ids,
					   unsigned int PDFfactor) const {
  Complex gvvh(0.,0.);
  getCouplings(gvvh,ids);
  double pre = norm(gvvh);
  if(sm_) {
    switch (PDFfactor) {
    case 0:
      return exp(2.*r/(pre));
    case 1:
    case 2:
    case 3:
    default:
      throw Exception() << "OneOneZeroEWSplitFn::invIntegOverP() invalid PDFfactor = "
  		      << PDFfactor << Exception::runerror;
    }
  }
  else {
    double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
    if(m0!=0.&&m1!=0.) {
      if(m0<m1) pre /= sqr(m0);
      else pre /= sqr(m1);
    }
    else if(m1!=0.) pre /= sqr(m1);
    else if(m0!=0.) pre /= sqr(m0);
    switch (PDFfactor) {
    case 0:
      return exp(2.*r/pre)/(1.+exp(2.*r/pre));
    case 1:
    case 2:
    case 3:
    default:
      throw Exception() << "OneOneZeroEWSplitFn::integOverP() invalid PDFfactor = "
  		      << PDFfactor << Exception::runerror;
    }
  }
}

bool OneOneZeroEWSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3)
    return false;
  if(_couplingValueIm==0&&_couplingValueRe==0) {
    if(abs(ids[0]->id())==ParticleID::Wplus && ids[2]->iSpin()==PDT::Spin0)
      return true;
    else if(ids[0]->id()==ParticleID::Z0 && ids[2]->iSpin()==PDT::Spin0)
      return true;
  }
  else { // BSM branching
    if(ids[0]->iSpin()==PDT::Spin1 && ids[1]->iSpin()==PDT::Spin1 && ids[2]->iSpin()==PDT::Spin0)
      return true;
  }
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
  Complex gvvh(0.,0.);
  getCouplings(gvvh,ids);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  Complex phase  = exp(Complex(0.,1.)*phi);
  Complex cphase = conj(phase);
  if(sm_) {
    double mVt = _theSM->mass(t,getParticleData(ids[0]->id()))/sqrt(t);
    double mHt = _theSM->mass(t,getParticleData(ids[2]->id()))/sqrt(t);
    double sqrtmass = sqrt(-(sqr(mVt)*sqr(1.-z))-sqr(mHt)*z+(1.-z)*z);
    // assign kernel
    (*kernal)(0,0,0) = -gvvh*mVt/sqrt(2.); 
    (*kernal)(0,1,0) = -gvvh*sqrtmass*cphase/2.; 
    (*kernal)(0,2,0) = 0.; 
    (*kernal)(1,0,0) = gvvh*sqrtmass*phase/2./z; 
    (*kernal)(1,1,0) = 0.;
    (*kernal)(1,2,0) = -gvvh*sqrtmass*cphase/2./z; 
    (*kernal)(2,0,0) = 0.;
    (*kernal)(2,1,0) = gvvh*sqrtmass*phase/2.; 
    (*kernal)(2,2,0) = -gvvh*mVt/sqrt(2.); 
  }
  else {
    double m0 = ids[0]->mass()/GeV, m1 = ids[1]->mass()/GeV;
    double m0t = ids[0]->mass()/sqrt(t), m1t = ids[1]->mass()/sqrt(t), m2t = ids[2]->mass()/sqrt(t);
    double sqrtmass = sqrt(z*(1.-z)*sqr(m0t)-(1.-z)*sqr(m1t)-z*sqr(m2t)+z*(1.-z));
    if(m0!=0.&&m1!=0.) {
      if(m0<m1) {
        gvvh /= m0;
        double mR = m0/m1;
        (*kernal)(0,0,0) = -gvvh*m0t/sqrt(2.); 
        (*kernal)(0,1,0) = -gvvh*mR*sqrtmass*cphase/2.; 
        (*kernal)(1,0,0) = gvvh*sqrtmass*phase/2./z; 
        (*kernal)(1,2,0) = -gvvh*sqrtmass*cphase/2./z; 
        (*kernal)(2,1,0) = gvvh*mR*sqrtmass*phase/2.; 
        (*kernal)(2,2,0) = -gvvh*m0t/sqrt(2.); 
      }
      else {
        gvvh /= m1;
        double mR = m1/m0;
        (*kernal)(0,0,0) = -gvvh*m1t/sqrt(2.); 
        (*kernal)(0,1,0) = -gvvh*sqrtmass*cphase/2.; 
        (*kernal)(1,0,0) = gvvh*mR*sqrtmass*phase/2./z; 
        (*kernal)(1,2,0) = -gvvh*mR*sqrtmass*cphase/2./z; 
        (*kernal)(2,1,0) = gvvh*sqrtmass*phase/2.; 
        (*kernal)(2,2,0) = -gvvh*m1t/sqrt(2.); 
      }
    }
    else if(m1!=0.) {
      gvvh /= m1;
      (*kernal)(0,0,0) = -gvvh*m1t/sqrt(2.); 
      (*kernal)(0,1,0) = -gvvh*sqrtmass*cphase/2.; 
      (*kernal)(1,0,0) = 0.;
      (*kernal)(1,2,0) = 0.;
      (*kernal)(2,1,0) = gvvh*sqrtmass*phase/2.; 
      (*kernal)(2,2,0) = -gvvh*m1t/sqrt(2.); 
    }
    else if(m0!=0.) {
      gvvh /= m0;
      (*kernal)(0,0,0) = -gvvh*m0t/sqrt(2.); 
      (*kernal)(0,1,0) = 0.;
      (*kernal)(1,0,0) = gvvh*sqrtmass*phase/2./z; 
      (*kernal)(1,2,0) = -gvvh*sqrtmass*cphase/2./z; 
      (*kernal)(2,1,0) = 0.;
      (*kernal)(2,2,0) = -gvvh*m0t/sqrt(2.); 
    }
    else {
      (*kernal)(0,0,0) = -gvvh*m0t/sqrt(2.); 
      (*kernal)(0,1,0) = 0.;
      (*kernal)(1,0,0) = 0.;
      (*kernal)(1,2,0) = 0.;
      (*kernal)(2,1,0) = 0.;
      (*kernal)(2,2,0) = -gvvh*m0t/sqrt(2.); 
    }
    (*kernal)(0,2,0) = 0.; 
    (*kernal)(1,1,0) = 0.;
    (*kernal)(2,0,0) = 0.;
  }
  // return the answer
  return kernal;
}
