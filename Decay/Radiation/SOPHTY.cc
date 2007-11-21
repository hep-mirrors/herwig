// -*- C++ -*-
//
// SOPHTY.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SOPHTY class.
//

#include "SOPHTY.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "FFDipole.h"
#include "IFDipole.h"

using namespace Herwig;

void SOPHTY::persistentOutput(PersistentOStream & os) const {
  os << _ffdipole << _ifdipole;
}

void SOPHTY::persistentInput(PersistentIStream & is, int) {
  is >> _ffdipole >> _ifdipole;
}

ClassDescription<SOPHTY> SOPHTY::initSOPHTY;
// Definition of the static class description member.

void SOPHTY::Init() {
  
  static ClassDocumentation<SOPHTY> documentation
    ("The SOPHTY class implements photon radiation in decays",
     "QED in particle decays was generated using the approach described in "
     "\\cite{Hamilton:2006xz}.",
     "\\bibitem{Hamilton:2006xz} K.~Hamilton and P.~Richardson,"
     "JHEP 07 (2006) 010.");
  
  static Reference<SOPHTY,FFDipole> interfaceFFDipole
    ("FFDipole",
     "The final-final dipole",
     &SOPHTY::_ffdipole, false, false, true, false, false);
  
  static Reference<SOPHTY,IFDipole> interfaceIFDipole
    ("IFDipole",
     "_ifdipole",
     &SOPHTY::_ifdipole, false, false, true, false, false);
}

ParticleVector SOPHTY::generatePhotons(const Particle & p,ParticleVector children) {
  if(children.size()!=2) return children;
  useMe();
  // final-final dipole
  if(p.dataPtr()->iCharge()==0) {
    if(children[0]->dataPtr()->iCharge()!=0&&
       children[1]->dataPtr()->iCharge()!=0)
      return _ffdipole->generatePhotons(p,children);
    else
      return children;
  }
  // initial final dipole
  else {
    if((children[0]->dataPtr()->iCharge()==0&&
	children[1]->dataPtr()->iCharge()!=0)||
       (children[0]->dataPtr()->iCharge()!=0&&
	children[1]->dataPtr()->iCharge()==0))
      return _ifdipole->generatePhotons(p,children);
    else
      return children;
  }
}
