// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PTVeto class.
//

#include "PTVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PTVeto.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PTVeto::~PTVeto() {}

void PTVeto::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << ounit(_maxPT,GeV) << ounit(_minPT,GeV) << _vetoTimelike << _vetoSpacelike << _vetoDecay;
}

void PTVeto::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> iunit(_maxPT,GeV) >> iunit(_minPT,GeV) >> _vetoTimelike >> _vetoSpacelike >> _vetoDecay;
}

ClassDescription<PTVeto> PTVeto::initPTVeto;
// Definition of the static class description member.

void PTVeto::Init() {

  
  static Parameter<PTVeto,Energy> interfacemaxPT
    ("maxPT",
     "Maximum allowed PT",
     &PTVeto::_maxPT, GeV, 100.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Parameter<PTVeto,Energy> interfaceminPT
    ("minPT",
     "Minimum allowed PT",
     &PTVeto::_minPT, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);



  static Switch<PTVeto,bool> interfacevetoTimelike
    ("vetoTimelike",
     "Veto timelike showering",
     &PTVeto::_vetoTimelike, true, false, false);
  static SwitchOption interfacevetoTimelikevetoTimelikeOn
    (interfacevetoTimelike,
     "vetoTimelikeOn",
     "Veto timelike showering",
     true);
  static SwitchOption interfacevetoTimelikevetoTimelikeOff
    (interfacevetoTimelike,
     "vetoTimelikeOff",
     "Do not veto timelike showering",
     false);

  static Switch<PTVeto,bool> interfacevetoSpacelike
    ("vetoSpacelike",
     "Veto spacelike showering",
     &PTVeto::_vetoSpacelike, true, false, false);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOn
    (interfacevetoSpacelike,
     "vetoSpacelikeOn",
     "Veto spacelike showering",
     true);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOff
    (interfacevetoSpacelike,
     "vetoSpacelikeOff",
     "Do not veto spacelike showering",
     false);

  static Switch<PTVeto,bool> interfacevetoDecay
    ("vetoDecay",
     "Veto decay showering",
     &PTVeto::_vetoDecay, true, false, false);
  static SwitchOption interfacevetoDecayvetoDecayOn
    (interfacevetoDecay,
     "vetoDecayOn",
     "Veto decay showering",
     true);
  static SwitchOption interfacevetoDecayvetoDecayOff
    (interfacevetoDecay,
     "vetoDecayOff",
     "Do not veto decay showering",
     false);


  static ClassDocumentation<PTVeto> documentation
    ("PT vetoing in shower");

}

