// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PomeronFlux class.
//

#include "PomeronFlux.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace ThePEG;
using namespace Herwig;

PomeronFlux::PomeronFlux() 
  : q2min_(ZERO)    , q2max_(10.*GeV2), 
    xiMin_(1.0e-7)  , xiMax_(1.0), 
    alfa0P_(1.104)  , alfapP_(0.06/GeV2), 
    betaP_(5.5/GeV2), normP_(1./GeV2),
    alfa0R_(0.5)    , alfapR_(0.3/GeV2), 
    betaR_(1.6/GeV2), normR_(1./GeV2),
    nR_(0.0013)     , PDFFit_(0)
{}

bool PomeronFlux::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == ParticleID::pplus );
}

cPDVector PomeronFlux::partons(tcPDPtr) const {
  cPDVector ret;
  ret.push_back(getParticleData( ParticleID::pomeron));
  ret.push_back(getParticleData( ParticleID::reggeon));
  return ret;
}

double PomeronFlux::xfl(tcPDPtr, tcPDPtr parton, Energy2 qq,
			double l, Energy2 ) const {
  // extra factor of Q^2 to cancel jacobian and sort out dimensions
  if(parton->id()==ParticleID::pomeron) 
    return normP_*qq*exp(-(betaP_ + 2.*alfapP_*l)*qq + 2.*(alfa0P_ - 1)*l);
  else if(parton->id()==ParticleID::reggeon)
    return nR_*normR_*qq*exp(-(betaR_ + 2.*alfapR_*l)*qq + 2.*(alfa0R_ - 1)*l);
  else {
    assert(false);
    return 0.;
  }
}

double PomeronFlux::xfvl(tcPDPtr, tcPDPtr, Energy2, double,
				   Energy2) const {
  // valence density is zero
  return 0.0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PomeronFlux,PDFBase>
describeHerwigPomeronFlux("Herwig::PomeronFlux", "HwPomeronFlux.so");

void PomeronFlux::Init() {

  static ClassDocumentation<PomeronFlux> documentation
    ("The PomeronFlux provides the PDF for a pomeron inside"
     " an incoming proton");

  static Parameter<PomeronFlux,Energy2> interfaceQ2Min
    ("Q2Min",
     "Minimum value of the magnitude of Q^2 for the pomeron/reggeon",
     &PomeronFlux::q2min_, GeV2, ZERO, ZERO, 100.0*GeV2,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,Energy2> interfaceQ2Max
    ("Q2Max",
     "Maximum value of the magnitude of Q^2 for the pomeron/reggeon",
     &PomeronFlux::q2max_, GeV2, 4.0*GeV2, ZERO, 100.0*GeV2,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,double> interfaceXiMin
    ("XiMin",
     "Minimum value of the pomeron/reggeon xi",
     &PomeronFlux::xiMin_,  1.0e-7 , 0., 1.,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,double> interfaceXiMax
    ("XiMax",
     "Maximum value of the pomeron/reggeon xi",
     &PomeronFlux::xiMax_, 1., 0., 1.,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,double> interfaceAlpha0Pomeron
    ("Alpha0Pomeron",
     "The coefficient of the intercept of alpha pomeron",
     &PomeronFlux::alfa0P_, 1.104, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<PomeronFlux,InvEnergy2> interfaceAlphaPrimePomeron
    ("AlphaPrimePomeron",
     "The coefficient of the t dependence of alpha of pomeron",
     &PomeronFlux::alfapP_, 1./GeV2, 0.06/GeV2, 0./GeV2, 1./GeV2,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,InvEnergy2> interfaceBPomeron
    ("BPomeron",
     "The coefficient of the t dependence of the exponential term of pomeron",
     &PomeronFlux::betaP_, 1./GeV2, 5.5/GeV2, 0./GeV2, 10./GeV2,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,double> interfaceAlpha0Reggeon
    ("Alpha0Reggeon",
     "The coefficient of the intercept of alpha reggeon",
     &PomeronFlux::alfa0R_, 0.5, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<PomeronFlux,InvEnergy2> interfaceAlphaPrimeReggeon
    ("AlphaPrimeReggeon",
     "The coefficient of the t dependence of alpha of reggeon",
     &PomeronFlux::alfapR_, 1./GeV2, 0.3/GeV2, 0./GeV2, 1./GeV2,
     false, false, Interface::limited);

  static Parameter<PomeronFlux,InvEnergy2> interfaceBReggeon
    ("BReggeon",
     "The coefficient of the t dependence of the exponential term reggeon",
     &PomeronFlux::betaR_, 1/GeV2, 1.6/GeV2, 0./GeV2, 10./GeV2,
     false, false, Interface::limited);
  
  static Parameter<PomeronFlux,double> interfacenR
    ("nR",
     "Reggeon flux factor",
     &PomeronFlux::nR_, 0.0017, 0.0, 0,
     true, false, Interface::lowerlim);

  static Switch<PomeronFlux,int> interfacePDFFit
    ("PDFFit",
     "The pomeron/regeon flux parameters are set according"
     " to choice of pomeron/regeon structure function fit. ",
     &PomeronFlux::PDFFit_, 0, false, false);
  static SwitchOption interfacePDFFitUser
    (interfacePDFFit,
     "User",
     "Default (set for fit Hera 2007) or user setting.",
      0);
  static SwitchOption interfacePDFFitPomeron2007
    (interfacePDFFit,
     "Pomeron2007",
     "Pomeron structure function fit is Hera 2007.",
      1);
  static SwitchOption interfacePDFFitPomeron2006A
    (interfacePDFFit,
     "Pomeron2006A",
     "Pomeron structure function fit is Hera 2006 A.",
     2);
  static SwitchOption interfacePDFFitPomeron2006B
    (interfacePDFFit,
     "Pomeron2006B",
     "Pomeron structure function fit is Hera 2006 B.",
     3);
}


double PomeronFlux::flattenScale(tcPDPtr proton, tcPDPtr parton, const PDFCuts & c,
	                         double l, double z, double & jacobian) const {
  double x = exp(-l);
  Energy2 qqmax = min(q2max_,0.25*sqr(x)*c.sMax());
  Energy2 qqmin = max(q2min_, sqr(proton->mass()*x)/(1-x));
  if(qqmin>=qqmax) {
    jacobian = 0.;
    return 0.;
  }
  InvEnergy2 k =ZERO;
  if(parton->id()==ParticleID::pomeron) 
    k = betaP_ + 2.*l*alfapP_; 
  else if(parton->id()==ParticleID::reggeon)
    k = betaR_ + 2.*l*alfapR_;
  else assert(false);
  double rho = exp(-k*qqmax) + z * ( exp(-k*qqmin) - exp(-k*qqmax) );
  Energy2 scale = -log(rho)/k;
  // jacobian factor (1/Q^2max_)(dQ^2/dz) / Q^2 (to make dimensionless)
  jacobian *=  (exp(-k*qqmin) - exp(-k*qqmax))/(rho*k*scale);
  // return Q^2/Q^2_max
  return scale/c.scaleMaxL(l);
}

double PomeronFlux::flattenL(tcPDPtr, tcPDPtr parton, const PDFCuts&,
			     double z, double & jacobian) const { 

  const double lMax = -log(xiMin_);
  const double lMin = -log(xiMax_);

  //  cout<<"lMin "<<lMin<<" lMax "<<lMax<<endl;

  if(parton->id()==ParticleID::pomeron) {
    jacobian *= lMax - lMin;
    return lMin + z*(lMax - lMin);
  }
  else if(parton->id()==ParticleID::reggeon) {
    double k = 2.*(alfa0R_ - 1.);
    double rho =exp(k*lMin) + z*(exp(k*lMax) - exp(k*lMin));
    jacobian *= (exp(k*lMax) - exp(k*lMin))/(rho*k);
    return log(rho)/k;
  }
  else {
    assert(false);
    return 0.;
  }
}

Energy2 PomeronFlux::intxFx(double x, Energy2 qqmin, Energy2 qqmax,
			    double alfa0, InvEnergy2 alfap,  InvEnergy2 beta) const {
  InvEnergy2 k =  beta - 2.*log(x)*alfap;
  return exp(-2.*(alfa0 - 1.)*log(x))*(exp(-qqmin*k) - exp(-qqmax*k))/k;
}

void PomeronFlux::persistentOutput(PersistentOStream & os) const {
  os << ounit(q2min_,GeV2) << ounit(q2max_,GeV2) << xiMin_ << xiMax_ << alfa0P_ 
     << ounit(alfapP_,1./GeV2) << ounit(betaP_,1./GeV2) << ounit(normP_,1./GeV2)
     << alfa0R_ << ounit(alfapR_,1./GeV2) << ounit(betaR_,1./GeV2) 
     << ounit(normR_,1./GeV2) << nR_ << PDFFit_;
}

void PomeronFlux::persistentInput(PersistentIStream & is, int) {
  is >> iunit(q2min_,GeV2) >> iunit(q2max_,GeV2) >> xiMin_ >> xiMax_ >> alfa0P_ 
     >> iunit(alfapP_,1./GeV2) >> iunit(betaP_,1./GeV2) >> iunit(normP_,1./GeV2)
     >> alfa0R_ >> iunit(alfapR_,1./GeV2) >> iunit(betaR_,1./GeV2) 
     >> iunit(normR_,1./GeV2) >> nR_ >> PDFFit_;
}

void PomeronFlux::doinit() {
  PDFBase::doinit();
  setFluxPar();
  // compute the normalisation factor
  double xp = 0.003;
  Energy2 tkinmin = sqr(getParticleData(ParticleID::pplus)->mass()*xp)/(1-xp);
  normP_ = 1./(intxFx(xp,tkinmin , 1.*GeV2, alfa0P_, alfapP_, betaP_ ));
  normR_ = 1./(intxFx(xp,tkinmin , 1.*GeV2, alfa0R_, alfapR_, betaR_ ));
}

void PomeronFlux::setFluxPar() {
  switch(PDFFit_){
  case 0:
    break;
  case 1:
    alfa0P_ = 1.104;
    nR_     = 0.0013;
    break;
  case 2:
    alfa0P_ = 1.118;
    nR_     = 0.0017;
    break;
  case 3:
    alfa0P_ = 1.111;
    nR_     = 0.0014;
    break;
  default:
    throw Exception() << "Invalid fit in PomeronFlux::setFluxPar()"
		      << Exception::runerror;
  }
}
