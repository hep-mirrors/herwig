// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GlobalParameters class.
//

#include "GlobalParameters.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>

using namespace Herwig;
// using namespace ThePEG;


GlobalParameters::~GlobalParameters() {}


void GlobalParameters::persistentOutput(PersistentOStream & os) const {
  os << _effectiveGluonMass 
     << _hadronizationScale
     << _stringFragmentationMode 
     << _softUnderlyingEventMode
     << _minVirtuality2
     << _maxDisplacement;
}


void GlobalParameters::persistentInput(PersistentIStream & is, int) {
  is >> _effectiveGluonMass 
     >> _hadronizationScale
     >> _stringFragmentationMode 
     >> _softUnderlyingEventMode
     >> _minVirtuality2
     >> _maxDisplacement;
}


ClassDescription<GlobalParameters> GlobalParameters::initGlobalParameters;
// Definition of the static class description member.


void GlobalParameters::Init() {

  static ClassDocumentation<GlobalParameters> documentation
    ("Class that stores Herwig global paramters");

  static Parameter<GlobalParameters,Energy> interfaceEffectiveGluonMass 
    ("EffectiveGluonMass",
     "Nonperturbative effective gluon mass  (unit [GeV])",
     &GlobalParameters::_effectiveGluonMass, 
     GeV, 0.750*GeV, 0.0*GeV, 10.0*GeV);

  static Parameter<GlobalParameters,Energy> interfaceHadronizationScale 
    ("HadronizationScale",
     "Hadronization scale, above which the decay affect the multi-scale showering  (unit [GeV])",
     &GlobalParameters::_hadronizationScale, 
     GeV, 0.5*GeV, 0.0*GeV, 10.0*GeV);

 static Switch<GlobalParameters, int> interfaceSoftUnderlyingEventMode
    ("OnOffSoftUnderlyingEventMode",
     "Choice of the soft underlying event switch mode.",
     &GlobalParameters::_softUnderlyingEventMode, 0, false, false);
  static SwitchOption interfaceSoftUnderlyingEventMode0
    (interfaceSoftUnderlyingEventMode,"SoftUnderlyingEvent-OFF", 
     "soft underlying event is OFF", 0);
  static SwitchOption interfaceSoftUnderlyingEventMode1
    (interfaceSoftUnderlyingEventMode,"SoftUnderlyingEvent-ON",
     "soft underlying event is ON", 1);

 static Switch<GlobalParameters, int> interfaceStringFragmentationMode
    ("OnOffThePEGStringFragmentationMode",
     "Choice of the ThePEG string fragmentation switch mode.",
     &GlobalParameters::_stringFragmentationMode, 0, false, false);
  static SwitchOption interfaceStringFragmentationMode0
    (interfaceStringFragmentationMode,"StringFragmentation-OFF", 
     "ThePEG string fragmentation is OFF", 0);
  static SwitchOption interfaceStringFragmentationMode1
    (interfaceStringFragmentationMode,"StringFragmentation-ON",
     "ThePEG string fragmentation is ON", 1);

  static Parameter<GlobalParameters,Energy2> interfaceMinVirtuality2 
    ("MinVirtuality2",
     "Minimum virtuality^2 of partons to use in calculating distances  (unit [GeV2]).",
     &GlobalParameters::_minVirtuality2, GeV2, 0.1*GeV2, 0.0*GeV2, 10.0*GeV2);

  static Parameter<GlobalParameters,Length> interfaceMaxDisplacement 
    ("MaxDisplacement",
     "Maximum displacement that is allowed for a particle  (unit [millimeter]).",
     &GlobalParameters::_maxDisplacement, millimeter, 1.0e-10*millimeter, 
     0.0*millimeter, 1.0e-9*millimeter);

}



