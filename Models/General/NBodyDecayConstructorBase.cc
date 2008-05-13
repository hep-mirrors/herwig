// -*- C++ -*-
//
// NBodyDecayConstructorBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NBodyDecayConstructorBase class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig; 
using namespace ThePEG;

void NBodyDecayConstructorBase::persistentOutput(PersistentOStream & os ) const {
  os << _init << _iteration << _points << _info;
}

void NBodyDecayConstructorBase::persistentInput(PersistentIStream & is , int) {
  is >> _init >> _iteration >> _points >> _info;
}

AbstractClassDescription<NBodyDecayConstructorBase> 
NBodyDecayConstructorBase::initNBodyDecayConstructorBase;
// Definition of the static class description member.

void NBodyDecayConstructorBase::Init() {

  static ClassDocumentation<NBodyDecayConstructorBase> documentation
    ("The NBodyDecayConstructorBase class is the base class for the automatic"
     "construction of the decay modes");
  
  static Switch<NBodyDecayConstructorBase,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &NBodyDecayConstructorBase::_init, true, false, false);
  static SwitchOption interfaceInitializeDecayersInitializeDecayersOn
    (interfaceInitializeDecayers,
     "Yes",
     "Initialize new decayers to find max weights",
     true);
  static SwitchOption interfaceInitializeDecayersoff
    (interfaceInitializeDecayers,
     "No",
     "Use supplied weights for integration",
     false);
  
  static Parameter<NBodyDecayConstructorBase,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &NBodyDecayConstructorBase::_iteration, 1, 0, 10,
     false, false, true);

  static Parameter<NBodyDecayConstructorBase,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &NBodyDecayConstructorBase::_points, 1000, 100, 100000000,
     false, false, true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceOutputInfo
    ("OutputInfo",
     "Whether to output information about the decayers",
     &NBodyDecayConstructorBase::_info, false, false, false);
  static SwitchOption interfaceOutputInfoOff
    (interfaceOutputInfo,
     "No",
     "Do not output information regarding the created decayers",
     false);
  static SwitchOption interfaceOutputInfoOn
    (interfaceOutputInfo,
     "Yes",
     "Output information regarding the decayers",
     true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceCreateDecayModes
    ("CreateDecayModes",
     "Whether to create the ThePEG::DecayMode objects as well as the decayers",
     &NBodyDecayConstructorBase::_createmodes, true, false, false);
  static SwitchOption interfaceCreateDecayModesOn
    (interfaceCreateDecayModes,
     "Yes",
     "Create the ThePEG::DecayMode objects",
     true);
  static SwitchOption interfaceCreateDecayModesOff
    (interfaceCreateDecayModes,
     "No",
     "Only create the Decayer objects",
     false);
}

void NBodyDecayConstructorBase::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  //Need width and branching ratios for all currently created decay modes
  PDPtr parent = const_ptr_cast<PDPtr>(dm->parent());
  DecaySet modes = parent->decayModes();
  if( modes.empty() ) return;
  double dmbrat(0.);
  if( modes.size() == 1 ) {
    parent->width(pwidth);
    parent->widthCut(5.*pwidth);
    dmbrat = 1.;
  }
  else {
    Energy currentwidth(parent->width());
    Energy newWidth(currentwidth + pwidth);
    parent->width(newWidth);
    parent->widthCut(5.*newWidth);
    //need to reweight current branching fractions if there are any
    for(DecaySet::const_iterator dit = modes.begin(); 
	dit != modes.end(); ++dit) {
      if( **dit == *dm || !(**dit).on() ) continue; 
      double newbrat = ((**dit).brat())*currentwidth/newWidth;
      ostringstream brf;
      brf << newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    //set brat for current mode
    dmbrat = pwidth/newWidth;
  }
  ostringstream br;
  br << dmbrat;
  generator()->preinitInterface(dm, "BranchingRatio",
				"set", br.str());
}

void NBodyDecayConstructorBase::setDecayerInterfaces(string fullname) const {
  if( initialize() ) {
    ostringstream value;
    value << initialize();
    generator()->preinitInterface(fullname, "Initialize", "set",
				  value.str());
    value.str("");
    value << iteration();
    generator()->preinitInterface(fullname, "Iteration", "set",
				  value.str());
    value.str("");
    value << points();
    generator()->preinitInterface(fullname, "Points", "set",
				  value.str());
  }
  string outputmodes;
  if( info() ) outputmodes = string("Output");
  else outputmodes = string("NoOutput");
  generator()->preinitInterface(fullname, "OutputModes", "set",
				outputmodes);
}
