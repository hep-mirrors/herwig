// -*- C++ -*-
//
// MPIPDF.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIPDF class.
//

#include "MPIPDF.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MPIPDF.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert> 

using namespace Herwig;
using namespace ThePEG;

MPIPDF::~MPIPDF() {}

bool MPIPDF::canHandleParticle(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->canHandleParticle(particle);
}

cPDVector MPIPDF::partons(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->partons(particle);
}

double MPIPDF::xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double l, Energy2 particleScale) const {
  // returns the density with removed valence part.
  assert(thePDF);
  return thePDF->xfl(particle, parton, partonScale, l, particleScale)-
    thePDF->xfvl(particle, parton, partonScale, l, particleScale);
}

double MPIPDF::xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale) const {
  // Here we should return the actual valence density.
  assert(thePDF);
  return thePDF->xfvl(particle, parton, partonScale, l, particleScale);
}


void MPIPDF::persistentOutput(PersistentOStream & os) const {
  os << thePDF;
}

void MPIPDF::persistentInput(PersistentIStream & is, int) {
  is >> thePDF;
}

ClassDescription<MPIPDF> MPIPDF::initMPIPDF;
// Definition of the static class description member.

void MPIPDF::Init() {

  static ClassDocumentation<MPIPDF> documentation
    ("There is no documentation for the MPIPDF class");

}

