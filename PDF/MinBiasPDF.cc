// -*- C++ -*-
//
// MinBiasPDF.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MinBiasPDF class.
//

#include "MinBiasPDF.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MinBiasPDF.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert> 

using namespace Herwig;
using namespace ThePEG;

bool MinBiasPDF::canHandleParticle(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->canHandleParticle(particle);
}

cPDVector MinBiasPDF::partons(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->partons(particle);
}

double MinBiasPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double eps, Energy2 particleScale) const {
  assert(thePDF);
  return thePDF->xfvx(particle, parton, partonScale, x, eps, particleScale);
}

double MinBiasPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps, Energy2 particleScale) const {
  assert(thePDF);
  return thePDF->xfvx(particle, parton, partonScale, x, eps, particleScale);
}


void MinBiasPDF::persistentOutput(PersistentOStream & os) const {
  os << thePDF;
}

void MinBiasPDF::persistentInput(PersistentIStream & is, int) {
  is >> thePDF;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MinBiasPDF,PDFBase>
describeHerwigMinBiasPDF("Herwig::MinBiasPDF", "HwShower.so");

void MinBiasPDF::Init() {

  static ClassDocumentation<MinBiasPDF> documentation
    ("MinBiasPDF is used to modify any given parton density. Currently it only returns the valence part");

  
  static Reference<MinBiasPDF,PDFBase> interfacePDF
    ("PDF",
     "pointer to the pdf, which will be modified",
     &MinBiasPDF::thePDF, false, false, true, false);
}

