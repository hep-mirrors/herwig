// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenHalfHalfOneSplitFn class.
//

#include "HiddenHalfHalfOneSplitFn.h"
#include "HiddenValleyModel.h"
#include "AdditionalGaugeParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerKinematics.h"

using namespace Herwig;

void HiddenHalfHalfOneSplitFn::doinit() {
  // get the model for parameters
  tcHiddenValleyPtr model = dynamic_ptr_cast<tcHiddenValleyPtr>
    (generator()->standardModel());
  if(!model) throw InitException() << "Must be using the HiddenValleyModel"
				   << " in HiddenHalfHalfOneSplitFn::doinit()" 
				   << Exception::runerror;
  assert(interactionType()!=ShowerInteraction::UNDEFINED);
  if(colourFactor()>0.) return;
  // compute the colour factors if need
  if(colourStructure()==TripletTripletOctet) {
    colourFactor(model->CF());
  }
  else if(colourStructure()==OctetOctetOctet) {
    colourFactor(model->CA());
  }
  else if(colourStructure()==OctetTripletTriplet) {
    colourFactor(model->TR());
  }
  else if(colourStructure()==TripletOctetTriplet) {
    colourFactor(model->CF());
  }
  else {
    assert(false);
  }
  HalfHalfOneSplitFn::doinit();
}


IBPtr HiddenHalfHalfOneSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr HiddenHalfHalfOneSplitFn::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<HiddenHalfHalfOneSplitFn> HiddenHalfHalfOneSplitFn::initHiddenHalfHalfOneSplitFn;
// Definition of the static class description member.

void HiddenHalfHalfOneSplitFn::Init() {

  static ClassDocumentation<HiddenHalfHalfOneSplitFn> documentation
    ("There is no documentation for the HiddenHalfHalfOneSplitFn class");

}

bool HiddenHalfHalfOneSplitFn::checkColours(const IdList & ids) const {
  tcAdditionalGaugeParticleDataPtr 
    pd[3]={dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(getParticleData(ids[0])),
	   dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(getParticleData(ids[1])),
	   dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(getParticleData(ids[2]))};
  if(colourStructure()==TripletTripletOctet) {
    if(ids[0]!=ids[1]) return false;
    if((pd[0]->hiddenColour()==HiddenPDT::HiddenColourFundamental||
	pd[0]->hiddenColour()==HiddenPDT::HiddenColourAntiFundamental) &&
       pd[2]->hiddenColour()==HiddenPDT::HiddenColourAdjoint) return true;
    return false;
  }
  else if(colourStructure()==OctetOctetOctet) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(pd[ix]->hiddenColour()!=HiddenPDT::HiddenColourAdjoint) return false;
    }
    return true;
  }
  else if(colourStructure()==OctetTripletTriplet) {
    if(pd[0]->hiddenColour()!=HiddenPDT::HiddenColourAdjoint) return false;
    if(pd[1]->hiddenColour()==HiddenPDT::HiddenColourFundamental&&
       pd[2]->hiddenColour()==HiddenPDT::HiddenColourAntiFundamental)
      return true;
    if(pd[1]->hiddenColour()==HiddenPDT::HiddenColourAntiFundamental&&
       pd[2]->hiddenColour()==HiddenPDT::HiddenColourFundamental)
      return true;
    return false;
  }
  else if(colourStructure()==TripletOctetTriplet) {
    if(ids[0]!=ids[2]) return false;
    if((pd[0]->hiddenColour()==HiddenPDT::HiddenColourFundamental||
	pd[0]->hiddenColour()==HiddenPDT::HiddenColourAntiFundamental) &&
       pd[1]->hiddenColour()==HiddenPDT::HiddenColourAdjoint) return true;
    return false;
  }
  else {
    assert(false);
  }
  return false;
}
