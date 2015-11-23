// -*- C++ -*-
//
// WeakDecayCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakDecayCurrent class.
//

#include "WeakDecayCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void WeakDecayCurrent::persistentOutput(PersistentOStream & os) const {
  os << _quark << _antiquark << _numbermodes;
}

void WeakDecayCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _quark >> _antiquark >> _numbermodes;
}

AbstractClassDescription<WeakDecayCurrent> WeakDecayCurrent::initWeakDecayCurrent;
// Definition of the static class description member.

void WeakDecayCurrent::Init() {

  static ClassDocumentation<WeakDecayCurrent> documentation
    ("The WeakDecayCurrent class is the basse class for the"
     " implementation of hadronic currents in weak decays.");

  static ParVector<WeakDecayCurrent,int> interfaceQuark
    ("Quark",
     "The PDG code for the quark.",
     &WeakDecayCurrent::_quark,
     0, 0, 0, 0, 16, false, false, true);

  static ParVector<WeakDecayCurrent,int> interfaceAntiQuark
    ("AntiQuark",
     "The PDG code for the antiquark.",
     &WeakDecayCurrent::_antiquark,
     0, 0, 0, -16, 0, false, false, true);
}

void WeakDecayCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) {
    output << "update decayers set parameters=\"";
  }
  if(create) {
    output << "create Herwig::WeakDecayCurrent " << name() << " \n";
  }
  for(unsigned int ix=0;ix<_quark.size();++ix) {
    if(ix<_numbermodes) {
      output << "newdef " << name() << ":Quark "     
	     << ix << "  " << _quark[ix]     << endl;
      output << "newdef " << name() << ":AntiQuark " 
	     << ix << "  " << _antiquark[ix] << endl;
    }
    else {
      output << "insert "  << name() << ":Quark "     
	     << ix << "  " << _quark[ix]     << endl;
      output << "insert "  << name() << ":AntiQuark " 
	     << ix << "  " << _antiquark[ix] << endl;
    }
  }
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
