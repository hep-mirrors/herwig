// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AdditionalGaugeParticleData class.
//

#include "AdditionalGaugeParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

AdditionalGaugeParticleData::AdditionalGaugeParticleData(long newId, string newPDGName)
  : ConstituentParticleData(newId, newPDGName), 
    hiddenColour_(HiddenPDT::HiddenColourNeutral)
{}

PDPtr AdditionalGaugeParticleData::
Create(long newId, string newPDGName) {
  return new_ptr(AdditionalGaugeParticleData(newId, newPDGName));
}

PDPair AdditionalGaugeParticleData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(AdditionalGaugeParticleData(newId, newPDGName));
  pap.second = new_ptr(AdditionalGaugeParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

void AdditionalGaugeParticleData::readSetup(istream & is) {
  ConstituentParticleData::readSetup(is);
  is >> ienum(hiddenColour_);
}

PDPtr AdditionalGaugeParticleData::pdclone() const {
  return new_ptr(*this);
}

void AdditionalGaugeParticleData::persistentOutput(PersistentOStream & os) const {
  os << oenum(hiddenColour_);
}

void AdditionalGaugeParticleData::persistentInput(PersistentIStream & is, int) {
  is >> ienum(hiddenColour_);
}

ClassDescription<AdditionalGaugeParticleData> 
AdditionalGaugeParticleData::initAdditionalGaugeParticleData;
// Definition of the static class description member.

void AdditionalGaugeParticleData::Init() {

  static ClassDocumentation<AdditionalGaugeParticleData> documentation
    ("The AdditionalGaugeParticleData class provides storage of the"
     " charge of the paritcle under the new guage group in hidden valley models");

  static Switch<AdditionalGaugeParticleData,HiddenPDT::HiddenColour> interfaceHiddenColour
    ("HiddenColour",
     "The charge under the new guage group",
     &AdditionalGaugeParticleData::hiddenColour_, HiddenPDT::HiddenColourNeutral,
     false, false);
  static SwitchOption interfaceHiddenColourNeutral
    (interfaceHiddenColour,
     "Neutral",
     "Uncharged under the new group",
     HiddenPDT::HiddenColourNeutral);
  static SwitchOption interfaceHiddenColourFundamental
    (interfaceHiddenColour,
     "Fundamental",
     "In the fundamental of the new group",
     HiddenPDT::HiddenColourFundamental);
  static SwitchOption interfaceHiddenColourAntiFundamental
    (interfaceHiddenColour,
     "AntiFundamental",
     "In the anti fundamental of the new gauge group",
     HiddenPDT::HiddenColourAntiFundamental);
  static SwitchOption interfaceHiddenColourAdjoint
    (interfaceHiddenColour,
     "Adjoint",
     "In the adjoint of the new gauge group",
     HiddenPDT::HiddenColourAdjoint);

}

