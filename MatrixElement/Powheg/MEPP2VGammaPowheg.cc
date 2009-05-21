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
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

MEPP2VGammaPowheg::MEPP2VGammaPowheg() 
  : _contrib(1), _scaleopt(0),
    _fixedScale(100.*GeV), _scaleFact(1.)
{}


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
  return MEPP2VGamma::dSigHatDR()*NLOweight();
}

double MEPP2VGammaPowheg::NLOweight() const {
  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
  return 1.;
}
