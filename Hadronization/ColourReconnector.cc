// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>

using namespace Herwig;

void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os << _clreco;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is >> _clreco;
}

ClassDescription<ColourReconnector> ColourReconnector::initColourReconnector;
// Definition of the static class description member.


void ColourReconnector::Init() {

  static ClassDocumentation<ColourReconnector> documentation
    ("This class is responsible of the colour reconnection.");


  static Switch<ColourReconnector,int> interfaceColourReconnection
    ("ColourReconnection",
     "Colour reconnections",
     &ColourReconnector::_clreco, 0, true, false);
  static SwitchOption interfaceColourReconnectionOff
    (interfaceColourReconnection,
     "Off",
     "Colour reconnections off",
     0);
  static SwitchOption interfaceColourReconnectionOn
    (interfaceColourReconnection,
     "On",
     "Colour reconnections on",
     1);
  
}


void ColourReconnector::rearrange(EventHandler &, 
				  ClusterVector &) 
   throw(Veto, Stop, Exception){
  if (_clreco != 0)
    throw Exception("Colour reconnection not implemented.",Exception::abortnow);
}
