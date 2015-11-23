// -*- C++ -*-
//
// SOPHTY.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SOPHTY class.
//

#include "SOPHTY.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "FFDipole.h"
#include "IFDipole.h"

using namespace Herwig;

void SOPHTY::persistentOutput(PersistentOStream & os) const {
  os << FFDipole_ << IFDipole_ << colouredOption_;
}

void SOPHTY::persistentInput(PersistentIStream & is, int) {
  is >> FFDipole_ >> IFDipole_ >> colouredOption_;
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
     &SOPHTY::FFDipole_, false, false, true, false, false);
  
  static Reference<SOPHTY,IFDipole> interfaceIFDipole
    ("IFDipole",
     "_ifdipole",
     &SOPHTY::IFDipole_, false, false, true, false, false);

  static Switch<SOPHTY,unsigned int> interfaceColouredTreatment
    ("ColouredTreatment",
     "Option for the treatment of QED radiation in decays involving coloured particles.",
     &SOPHTY::colouredOption_, 0, false, false);
  static SwitchOption interfaceColouredTreatmentNone
    (interfaceColouredTreatment,
     "None",
     "Generate no QED radiation to avoid problems with the interplay"
     " of QCD and QED radiation",
     0);
  static SwitchOption interfaceColouredTreatmentRadiation
    (interfaceColouredTreatment,
     "Radiation",
     "Generate radiation from the coloured particles.",
     1);

}

ParticleVector SOPHTY::generatePhotons(const Particle & p,ParticleVector children,
				       tDecayIntegratorPtr decayer) {
  if ( children.size() != 2 ) return children;
  // if not generating radiation from coloured particles
  // return if there are any coloured particles
  if(colouredOption_==0) {
    bool coloured = p.dataPtr()->coloured();
    for(unsigned int ix=0;ix<children.size();++ix) {
      coloured |= children[ix]->dataPtr()->coloured();
    }
    if(coloured) return children;
  }
  useMe();
  // final-final dipole
  if(p.dataPtr()->iCharge()==0) {
    if(children[0]->dataPtr()->iCharge()!=0&&
       children[1]->dataPtr()->iCharge()!=0)
      return FFDipole_->generatePhotons(p,children,decayer);
    else
      return children;
  }
  // initial final dipole
  else {
    if((children[0]->dataPtr()->iCharge()==0&&
	children[1]->dataPtr()->iCharge()!=0)||
       (children[0]->dataPtr()->iCharge()!=0&&
	children[1]->dataPtr()->iCharge()==0))
      return IFDipole_->generatePhotons(p,children);
    else
      return children;
  }
}
