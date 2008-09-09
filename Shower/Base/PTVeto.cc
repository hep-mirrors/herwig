// -*- C++ -*-
//
// PTVeto.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PTVeto class.
//

#include "PTVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PTVeto::persistentOutput(PersistentOStream & os) const {
  os << ounit(_maxPT,GeV) << ounit(_minPT,GeV) << _vetoTimelike << _vetoSpacelike;
}

void PTVeto::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_maxPT,GeV) >> iunit(_minPT,GeV) >> _vetoTimelike >> _vetoSpacelike;
}

ClassDescription<PTVeto> PTVeto::initPTVeto;
// Definition of the static class description member.

void PTVeto::Init() {

  
  static Parameter<PTVeto,Energy> interfacemaxPT
    ("MaxPT",
     "Maximum allowed PT",
     &PTVeto::_maxPT, GeV, 100.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Parameter<PTVeto,Energy> interfaceminPT
    ("MinPT",
     "Minimum allowed PT",
     &PTVeto::_minPT, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);



  static Switch<PTVeto,bool> interfacevetoTimelike
    ("VetoTimelike",
     "Veto timelike showering",
     &PTVeto::_vetoTimelike, true, false, false);
  static SwitchOption interfacevetoTimelikevetoTimelikeOn
    (interfacevetoTimelike,
     "Yes",
     "Veto timelike showering",
     true);
  static SwitchOption interfacevetoTimelikevetoTimelikeOff
    (interfacevetoTimelike,
     "No",
     "Do not veto timelike showering",
     false);

  static Switch<PTVeto,bool> interfacevetoSpacelike
    ("VetoSpacelike",
     "Veto spacelike showering",
     &PTVeto::_vetoSpacelike, true, false, false);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOn
    (interfacevetoSpacelike,
     "Yes",
     "Veto spacelike showering",
     true);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOff
    (interfacevetoSpacelike,
     "No",
     "Do not veto spacelike showering",
     false);


  static ClassDocumentation<PTVeto> documentation
    ("PT vetoing in shower");

}

