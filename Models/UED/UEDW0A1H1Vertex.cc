// -*- C++ -*-
//
// UEDW0A1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0A1H1Vertex class.
//

#include "UEDW0A1H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

void UEDW0A1H1Vertex::doinit() {
  addToList( 24, 5100036, -5100037);
  addToList(-24, 5100036,  5100037);
  VSSVertex::doinit();
  tUEDBasePtr UEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase) throw InitException() 
    << "UEDW0A1H1Vertex::doinit() - The pointer to "
    << "the UEDBase object is null!" << Exception::runerror;
  theMw2 = sqr(getParticleData(24)->mass());
  theMz2 = sqr(getParticleData(23)->mass());
  theR2 = sqr(UEDBase->compactRadius());
}

UEDW0A1H1Vertex::UEDW0A1H1Vertex() : theMw2(), theMz2(), theR2(), 
				     theq2Last(), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDW0A1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMw2,GeV2) << ounit(theMz2,GeV2) << ounit(theR2,1/GeV2);  
}

void UEDW0A1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMw2,GeV2) >> iunit(theMz2,GeV2) >> iunit(theR2,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDW0A1H1Vertex,VSSVertex>
describeHerwigUEDW0A1H1Vertex("Herwig::UEDW0A1H1Vertex", "HwUED.so");

void UEDW0A1H1Vertex::Init() {

  static ClassDocumentation<UEDW0A1H1Vertex> documentation
    ("The coupling of a SM W boson to a level-1 charged higgs and the "
     "level-1 heavy neutral higgs");

}

void UEDW0A1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long chiggs(0);
  if(abs(part1->id()) == ParticleID::Wplus)
    chiggs = (abs(part2->id()) == 5100037) ? part2->id() : part3->id();
  else if(abs(part2->id()) == ParticleID::Wplus)
    chiggs = (abs(part1->id()) == 5100037) ? part1->id() : part3->id();
  else if(abs(part3->id()) == ParticleID::Wplus)
    chiggs = (abs(part1->id()) == 5100037) ? part1->id() : part2->id();
  else {
    throw HelicityLogicalError() << "UEDW0A1H1Vertex::setCoupling - "
				 << "There is no SM W boson in this vertex"
				 << Exception::warning;
    return;
  }
  if(abs(chiggs) == 5100037) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = weakCoupling(q2);
      double mwRs = theMw2*theR2;
      double denom = sqrt( (1 + mwRs)*(1. + theMw2*theR2) );
      theCoupLast *= ( 0.5 + mwRs )/denom;
    }
    if(chiggs > 0) theCoupLast *= -1.;
    norm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDW0A1H1Vertex::setCoupling - "
				 << "There is an unknown particle in this " 
				 << "vertex " << chiggs
				 << Exception::runerror;
}

