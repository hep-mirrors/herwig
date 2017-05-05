// -*- C++ -*-
//
// SatPDF.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SatPDF class.
//

#include "SatPDF.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SatPDF.tcc"
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

SatPDF::~SatPDF() {}

bool SatPDF::canHandleParticle(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->canHandleParticle(particle);
}

cPDVector SatPDF::partons(tcPDPtr particle) const {
  assert(thePDF);
  return thePDF->partons(particle);
}

double SatPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double eps, Energy2 particleScale) const {
  assert(thePDF);
  double xUsed(x);
  if( x < theX0 )
    xUsed = theX0;

  double pdf(thePDF->xfx(particle, parton, partonScale, 
			 xUsed, eps, particleScale));

  if( x < theX0 )
    pdf *= pow(x/theX0, theExp);

  return pdf;
}

double SatPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps, Energy2 particleScale) const {
  // Here we should return the actual valence density.
  assert(thePDF);
  double xUsed(x);
  if( x < theX0 )
    xUsed = theX0;

  double pdf(thePDF->xfvx(particle, parton, partonScale, 
			 xUsed, eps, particleScale));
  if( x < theX0 )
    pdf *= pow(x/theX0, theExp);

  return pdf;
}


void SatPDF::persistentOutput(PersistentOStream & os) const {
  os << thePDF << theX0 << theExp;
}

void SatPDF::persistentInput(PersistentIStream & is, int) {
  is >> thePDF >> theX0 >> theExp;
}

ClassDescription<SatPDF> SatPDF::initSatPDF;
// Definition of the static class description member.

void SatPDF::Init() {

  static ClassDocumentation<SatPDF> documentation
    ("SatPDF is used to modify any given parton density for small values of x. That can be used to mimic saturation effects.");

  
  static Reference<SatPDF,PDFBase> interfacePDF
    ("PDF",
     "pointer to the pdf, which will be modified",
     &SatPDF::thePDF, false, false, true, false);

  
  static Parameter<SatPDF,double> interfaceX0
    ("X0",
     "For x values below this one the parametrisation: f(x) = f(x0) * (x/x0)**exp is used",
     &SatPDF::theX0, 1.E-4, 0.0, 1.0,
     true, false, Interface::limited);

  static Parameter<SatPDF,double> interfaceExp
    ("Exp",
     "For x values below X0 the parametrisation: f(x) = f(x0) * (x/x0)**exp is used",
     &SatPDF::theExp, 0.0, -10.0, 10.0,
     true, false, Interface::limited);


}

