// -*- C++ -*-
//
// BtoSGammaHadronicMass.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaHadronicMass class.
//

#include "BtoSGammaHadronicMass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void BtoSGammaHadronicMass::persistentOutput(PersistentOStream & os) const {
  os << ounit(_minMass,GeV) << ounit(_maxMass,GeV);
}

void BtoSGammaHadronicMass::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_minMass,GeV) >> iunit(_maxMass,GeV);
}

AbstractClassDescription<BtoSGammaHadronicMass> 
BtoSGammaHadronicMass::initBtoSGammaHadronicMass;
// Definition of the static class description member.

void BtoSGammaHadronicMass::Init() {

  static ClassDocumentation<BtoSGammaHadronicMass> documentation
    ("The BtoSGammaHadronicMass class is the base class for the implementation"
     " of models of the hadronic spectrum in B to s gamma decays.");

  static Parameter<BtoSGammaHadronicMass,Energy> interfaceMinimumMass
    ("MinimumMass",
     "The minimum value of the hadronic mass",
     &BtoSGammaHadronicMass::_minMass, GeV, 0.825*GeV, 0.825*GeV, 5.300*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaHadronicMass,Energy> interfaceMaximumMass
    ("MaximumMass",
     "The maximum value of the hadronic mass",
     &BtoSGammaHadronicMass::_maxMass, GeV, 5.300*GeV, 0.825*GeV, 5.300*GeV,
     false, false, Interface::limited);

}

void BtoSGammaHadronicMass::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BtoSGammaHadronicMass " 
		    << name() << " \n";
  output << "newdef " << name() << ":MinimumMass " << _minMass/GeV << " \n";
  output << "newdef " << name() << ":MaximumMass " << _maxMass/GeV << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
