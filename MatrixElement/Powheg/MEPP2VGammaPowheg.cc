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
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Utilities/Maths.h"


using namespace Herwig;
using Herwig::Math::ReLi2;

MEPP2VGammaPowheg::MEPP2VGammaPowheg() 
  :  _contrib(1), _scaleopt(0), _fixedScale(100.*GeV), _scaleFact(1.)
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
  _xa=  lastX1();
  _xb=  lastX2();
  // radative variables
  _x = r[nDim()];
  _z = r[nDim()+1];
  _phi = Constants::twopi*r[nDim()+2];
  return MEPP2VGamma::generateKinematics(r);
}

CrossSection MEPP2VGammaPowheg::dSigHatDR() const {
  return MEPP2VGamma::dSigHatDR()*NLOweight();
}

// definitions of H, F^V:
double MEPP2VGammaPowheg::Hfunc(Energy2 t, Energy2 s, Energy2 m2) const{  
  // this is eqn 16 of PRD 47, 940
  return sqr(Constants::pi) - sqr(log(s/m2)) + sqr(log(-t/s)) - sqr(log(-t/m2))
    -2.*Math::ReLi2(1.-s/m2)-2.*Math::ReLi2(1.-t/m2);
} 


double MEPP2VGammaPowheg::FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
  double y1 = 4. * ( 2.*sqr(s)/(t*u) + 2.*s/u + t/u ) * Hfunc(u,s,m2);
  double y2 = -8./3.*sqr(Constants::pi)*s*( 2.*s/t + t/s - u/s )/(t+u);
  double y3 = 4. * (6. - 10.*u/(t+u) - 10.*sqr(s)/(u*(t+u))
		    - 11.*s/u - 5.*t/u + 2.*s/(t+u) + s/(s+t));
  double y4 = -4.*log(s/m2) * (3.*t/u + 2.*s/u + 4.*s*(t+s)/(u*(t+u)) 
			       + 2.*t/u*sqr(s/(t+u)));
  double y5 = 4. * log(-u/m2)*((4.*s+u)/(s+t) + s*u/sqr(s+t));
  return y1 + y2 + y3 + y4 + y5;
}

double MEPP2VGammaPowheg::FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  // this is eqn 15 of PRD 47, 940, check by PR 23/6/09
  double y1 = 4. * ( 2.*sqr(s)/(t*u) + 2.*s/u + t/u ) * Hfunc(u,s,m2);
  double y2 = -8./3. * sqr(Constants::pi) * sqr(s)/(t*u);
  double y3 = 4. * ( 1. - 5.*sqr(s)/(t*u) - 11.*s/u-5.*t/u + 2.*s/(t+u) + s/(s+t) );
  double y4 = 4. * log( s/m2) * ( 4.*s/(t+u) + 2.*sqr(s/(t+u)) - 3.*sqr(s+t)/(t*u) );
  double y5 = 4. * log(-u/m2)*( (4.*s+u)/(s+t) + s*u/sqr(s+t) );
  return y1 + y2 + y3 + y4 + y5;
}
 
double MEPP2VGammaPowheg::FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  return mePartonData()[2]->id()==ParticleID::Z0 
    ? FZfunc(t,u,s,m2) : FWfunc(t,u,s,m2);
}

// ratio of NLO/LO 
double MEPP2VGammaPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
  // mass squared of the vector boson
  Energy2 MV2 = meMomenta()[2].m2();
  // charges of the quarks
  double charge0 = mePartonData()[0]->iCharge()/3.;
  double charge1 = mePartonData()[1]->iCharge()/3.;
  // the strong coupling
  double alphas= SM().alphaS(scale());
  // the different pieces
  // LO prefactor
  Energy4 smborn0 =   (sHat()*MV2 -tHat()*uHat() + 0.5*sqr(tHat()+uHat()));
  // LO O(e) piece (N.B. changed sign here PR 6/23/09
  Energy4 smborn1 = - (sHat()*MV2 -tHat()*uHat() +     sqr(tHat()+uHat()));
  // LO O(e^2) piece
  Energy4 smborn2 = 0.5*sqr(tHat()+uHat());
  // the piece of the virtual correction left after cancelling the poles
  double partep= 2. * (5. - sqr(Constants::pi)/3.0 ) 
    + 3.0 * smborn1/smborn0 + 2.0 * smborn2/smborn0;
  // the finite piece of the virtual correction
  // prefactors from the LO piece
  InvEnergy4 smbfact = 4.0 * sqr(charge0*tHat() + charge1*uHat())/(tHat()*uHat()*sqr(tHat()+uHat()));
  // the finite piece of the 1 loop with some prefactors removed
  double smloop = 0.5 * (charge0*tHat()+ charge1*uHat())/(tHat()+uHat())*
    (charge0*FVfunc(tHat(),uHat(),sHat(),MV2) + charge1*FVfunc(uHat(),tHat(),sHat(),MV2));
  // full virtual piece
  double partloop = smloop/(smborn0*smbfact);
  // alphaS prefactor
  double alfsfact = alphas*_CF/(2.0*Constants::pi);
  // virtual correction
  return 1. + alfsfact * (partloop + partep); 
//   return 1.;
}
