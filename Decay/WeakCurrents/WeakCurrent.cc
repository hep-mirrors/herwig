// -*- C++ -*-
//
// WeakCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakCurrent class.
//

#include "WeakCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG;

void WeakCurrent::persistentOutput(PersistentOStream & os) const {
  os << quark_ << antiquark_ << numberModes_;
}

void WeakCurrent::persistentInput(PersistentIStream & is, int) {
  is >> quark_ >> antiquark_ >> numberModes_;
}

void WeakCurrent::constructSpinInfo(ParticleVector decay) const {
  using namespace Helicity;
  for(unsigned int ix=0;ix<decay.size();++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<WeakCurrent,Interfaced>
describeHerwigWeakCurrent("Herwig::WeakCurrent", "Herwig.so");

void WeakCurrent::Init() {

  static ClassDocumentation<WeakCurrent> documentation
    ("The WeakCurrent class is the basse class for the"
     " implementation of hadronic currents in weak decays.");

  static ParVector<WeakCurrent,int> interfaceQuark
    ("Quark",
     "The PDG code for the quark.",
     &WeakCurrent::quark_,
     0, 0, 0, 0, 16, false, false, true);

  static ParVector<WeakCurrent,int> interfaceAntiQuark
    ("AntiQuark",
     "The PDG code for the antiquark.",
     &WeakCurrent::antiquark_,
     0, 0, 0, -16, 0, false, false, true);
}

void WeakCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) {
    output << "update decayers set parameters=\"";
  }
  if(create) {
    output << "create Herwig::WeakCurrent " << name() << " \n";
  }
  for(unsigned int ix=0;ix<quark_.size();++ix) {
    if(ix<numberModes_) {
      output << "newdef " << name() << ":Quark "     
	     << ix << "  " << quark_[ix]     << endl;
      output << "newdef " << name() << ":AntiQuark " 
	     << ix << "  " << antiquark_[ix] << endl;
    }
    else {
      output << "insert "  << name() << ":Quark "     
	     << ix << "  " << quark_[ix]     << endl;
      output << "insert "  << name() << ":AntiQuark " 
	     << ix << "  " << antiquark_[ix] << endl;
    }
  }
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
