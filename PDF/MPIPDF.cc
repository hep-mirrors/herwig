// -*- C++ -*-
//
// MPIPDF.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIPDF class.
//

#include "MPIPDF.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert> 

using namespace Herwig;
using namespace ThePEG;

IBPtr MPIPDF::clone() const {
  return new_ptr(*this);
}

IBPtr MPIPDF::fullclone() const {
  return new_ptr(*this);
}


bool MPIPDF::canHandleParticle(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->canHandleParticle(particle);
}

cPDVector MPIPDF::partons(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->partons(particle);
}

double MPIPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		   double x, double eps, Energy2 particleScale) const {
  // returns the density with removed valence part.
  assert(thePDF);
  return thePDF->xfsx(particle, parton, partonScale, 
		      x, eps, particleScale);
}

double MPIPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		    double x, double eps, Energy2 particleScale) const {
  // Here we should return the actual valence density.
  assert(thePDF);
  return thePDF->xfvx(particle, parton, partonScale, 
		      x, eps, particleScale);
}


void MPIPDF::persistentOutput(PersistentOStream & os) const {
  os << thePDF;
}

void MPIPDF::persistentInput(PersistentIStream & is, int) {
  is >> thePDF;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MPIPDF,PDFBase>
describeHerwigMPIPDF("Herwig::MPIPDF", "HwShower.so");

void MPIPDF::Init() {

  static ClassDocumentation<MPIPDF> documentation
    ("The MPIPDF class wraps other PDF classes to provide "
     "sea-quark-only PDFs.");

}

