// -*- C++ -*-
//
// HepMCFile.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HepMCFile class.
//

#include "HepMCFile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <HepMCHelper.h>

using namespace ThePEG;
using namespace Herwig;

void HepMCFile::analyze(tEventPtr event, long, int, int) {
  if (event->number() > _eventNumber) return;
  HepMC::GenEvent * hepmc = HepMCConverter<HepMC::GenEvent>::convert(*event);
  hepmc->print(_hepmcfile);
  delete hepmc;
}

void HepMCFile::persistentOutput(PersistentOStream & os) const {
  os << _eventNumber;
}

void HepMCFile::persistentInput(PersistentIStream & is, int) {
  is >> _eventNumber;
}


ClassDescription<HepMCFile> HepMCFile::initHepMCFile;
// Definition of the static class description member.

void HepMCFile::Init() {

  static ClassDocumentation<HepMCFile> documentation
    ("This analysis handler will output the event record in HepMC format.");

  static Parameter<HepMCFile,long> interfacePrintEvent
    ("PrintEvent",
     "The number of events that should be printed.",
     &HepMCFile::_eventNumber, 1, 1, 1,
     false, false, Interface::lowerlim);
}

