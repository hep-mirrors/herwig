// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VGammaPowheg class.
//

#include "MEPP2VGammaPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Utilities/GSLIntegrator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Utilities/Maths.h"


using namespace Herwig;
using Herwig::Math::ReLi2;

MEPP2VGammaPowheg::MEPP2VGammaPowheg() 
  :  _contrib(1), _scaleopt(0),
     _fixedScale(100.*GeV), _scaleFact(1.)
{}

void MEPP2VGammaPowheg::doinit() {
  // gluon ParticleData object
  //_gluon = getParticleData(ParticleID::g);
  // colour factors
  _CF = 4./3.; 
  _TR = 0.5;

  MEPP2VGamma::doinit();
}


IBPtr MEPP2VGammaPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2VGammaPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2VGammaPowheg::persistentOutput(PersistentOStream & os) const {
  os << _contrib << _scaleopt << ounit(_fixedScale,GeV) << _scaleFact;
}

void MEPP2VGammaPowheg::persistentInput(PersistentIStream & is, int) {
  is >> _contrib >> _scaleopt >> iunit(_fixedScale,GeV) >> _scaleFact;
}

ClassDescription<MEPP2VGammaPowheg> MEPP2VGammaPowheg::initMEPP2VGammaPowheg;
// Definition of the static class description member.

void MEPP2VGammaPowheg::Init() {

  static ClassDocumentation<MEPP2VGammaPowheg> documentation
    ("The MEPP2VGammaPowheg class implements the NLO matrix"
     " elements for q qbar -> W/Z gamma in the POWHEG scheme.");

   static Switch<MEPP2VGammaPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2VGammaPowheg::_contrib, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Switch<MEPP2VGammaPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the scale to be used",
     &MEPP2VGammaPowheg::_scaleopt, 0, false, false);
  static SwitchOption interfaceScaleOptionFixed
    (interfaceFactorizationScaleOption,
     "Fixed",
     "Use a fixed scale",
     0);
  static SwitchOption interfaceScaleOptionsHat
    (interfaceFactorizationScaleOption,
     "Dynamic",
     "Used sHat as the scale",
     1);

  static Parameter<MEPP2VGammaPowheg,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "The fixed scale to use if required",
     &MEPP2VGammaPowheg::_fixedScale, GeV, 100.0*GeV, 10.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEPP2VGammaPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2VGammaPowheg::_scaleFact, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
}

Energy2 MEPP2VGammaPowheg::scale() const {
  return _scaleopt == 0 ? sqr(_fixedScale) : _scaleFact*sHat();
}

int MEPP2VGammaPowheg::nDim() const {
  return MEPP2VGamma::nDim()+3;
}

bool MEPP2VGammaPowheg::generateKinematics(const double * r) {
  return MEPP2VGamma::generateKinematics(r);
}

CrossSection MEPP2VGammaPowheg::dSigHatDR() const {
  if(mePartonData()[3]->id()==ParticleID::gamma) {
    _idboson= 2;
    _p_photon= meMomenta()[3];
    _p_boson=  meMomenta()[2];}
  else{
    _idboson= 3;
    _p_photon= meMomenta()[2];
    _p_boson=  meMomenta()[3];
  }
  _xa=  lastX1();
  _xb=  lastX2();
  return MEPP2VGamma::dSigHatDR()*NLOweight();
}

// definition of Li_2 function:
double MEPP2VGammaPowheg::Li2p(double z){
  Li2Integrand integrand(z);
  return _integrator.value(integrand, 0.0, 1.0);
}


// definitions of H, F^V:
double MEPP2VGammaPowheg::Hfunc(Energy2 t, Energy2 s, Energy2 m2) const{  
  return sqr(Constants::pi)-sqr(log(s/m2))+sqr(log(-t/s))-sqr(log(-t/m2))
    -2.0*ReLi2(1.0-s/m2)-2.0*ReLi2(1.0-t/m2);
} 


double MEPP2VGammaPowheg::FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(Constants::pi)*s*(2.0*s/t +t/s -u/s)/(t+u);
    y3= 4.0*(6.0-10.0*u/(t+u)-10.0*s*s/(u*(t+u))-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= -4.0*log(s/m2)*(3.0*t/u+2.0*s/u+4.0*s*(t+s)/(u*(t+u))+2.0*t/u*sqr(s/(t+u)));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(Constants::pi)*s*s/(t*u);
    y3= 4.0*(1.0-5.0*s*s/(t*u)-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= 4.0*log(s/m2)*(4.0*s/(t+u)+2.0*sqr(s/(t+u))-3.0*sqr(s+t)/(t*u));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  if(mePartonData()[_idboson]->id()==ParticleID::Z0) {
    return FZfunc(t,u,s,m2);
  }
  else {
    return FWfunc(t,u,s,m2);
  }
}

// ratio of NLO/LO 
double MEPP2VGammaPowheg::NLOweight() const {
  //  double Pi(3.1415926);
  double partep, partloop, alfsfact;
  double smborn0, smborn1, smborn2, smbfact, smloop;
  Energy2 MV2;
  double charge0,charge1,alphae,sin2w,sh;
  //  double FV, FW, FZ;
  // If only leading order is required return 1:
   if(_contrib==0) return 1.;

   
  _ss= sHat();
  _tt= ( meMomenta()[0]- _p_boson).m2();
  _uu= ( meMomenta()[0]- _p_photon).m2();
  MV2= _p_boson.m2();
  //MV2= sqr(81.0*GeV);
  charge0= mePartonData()[0]->iCharge()/3.;
  charge1= mePartonData()[1]->iCharge()/3.;

  _alphas= SM().alphaS(scale());
  //_alphas=0.112;
  //_CF= 4.0/3.0;
  smbfact= 4.0*sqr(charge0*_tt + charge1*_uu)/(_tt*_uu*sqr(_tt+_uu))*GeV*GeV*GeV*GeV;
  smborn0= (_ss*MV2 -_tt*_uu +sqr(_tt+_uu)/2.0)/(GeV*GeV*GeV*GeV);
  smborn1= (_ss*MV2 -_tt*_uu +sqr(_tt+_uu))/(GeV*GeV*GeV*GeV);
  smborn2= (sqr(_tt+_uu)/2.0)/(GeV*GeV*GeV*GeV);
  partep= 2.0*(5.0-sqr(Constants::pi)/3.0)*1.0 +3.0*smborn1/smborn0 +2.0*smborn2/smborn0;

  smloop= (charge0*_tt+ charge1*_uu)*
       (charge0*FVfunc(_tt,_uu,_ss,MV2) +charge1*FVfunc(_uu,_tt,_ss,MV2))/(_tt+_uu)/2.0;
  partloop= smloop/(smborn0*smbfact);  
  alfsfact= _alphas*_CF/(2.0*Constants::pi);
  sh=_tt/(GeV*GeV);
  //cerr << "CF=" << _CF << " charge1=" << charge1 <<  " tHat=" << sh << "*GeV  virtual loop & dipole ratio = " 
  //     << alfsfact*(partep+partloop) << " virtual loop ratio = " <<  alfsfact*partloop << "\n";
  //return 1.+ alfsfact*(partep+partloop) ;
  //return 1.+_alphas*4.0/3.0/(2.0*Constants::pi)*(partep) ; 
  return 1.;
}
