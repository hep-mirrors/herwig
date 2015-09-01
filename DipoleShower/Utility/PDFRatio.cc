// -*- C++ -*-
//
// PDFRatio.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PDFRatio class.
//

#include "PDFRatio.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/PDF/HwRemDecayer.h"

using namespace Herwig;

PDFRatio::PDFRatio() 
  : HandlerBase(), 
    theValenceExtrapolation(0.7), theSeaExtrapolation(0.6),
    theFreezingScale(1.0*GeV) {}

PDFRatio::~PDFRatio() {}

IBPtr PDFRatio::clone() const {
  return new_ptr(*this);
}

IBPtr PDFRatio::fullclone() const {
  return new_ptr(*this);
}

double PDFRatio::operator() (const PDF& pdf,
			     Energy2 scale,
			     tcPDPtr from, tcPDPtr to,
			     double x, double z) const {

  if ( x/z > 1.0 )
    return 0.0;

  if ( scale < sqr(theFreezingScale) )
    scale = sqr(theFreezingScale);

  tcHwRemDecPtr remDec = 
    ShowerHandler::currentHandler()->remnantDecayer();

  const HwRemDecayer::HadronContent* theContent = 0;

  if ( remDec->content().first.hadron == pdf.particle() ) {
    theContent = &(remDec->content().first);
  } 

  if ( remDec->content().second.hadron == pdf.particle() ) {
    theContent = &(remDec->content().second);
  }

  assert(theContent);

  double exFrom =
    theContent->isValenceQuarkData(from) ?
    theValenceExtrapolation : theSeaExtrapolation;

  double exTo =
    theContent->isValenceQuarkData(to) ?
    theValenceExtrapolation : theSeaExtrapolation;

  double xTo = x/z;

  double fromPDF =
    x < exFrom ?
    pdf.xfx(from,scale,x) :
    ((1.-x)/(1.-exFrom)) * pdf.xfx(from,scale,exFrom);

  if ( fromPDF < 1e-8 )
    fromPDF = 0.0;

  double toPDF =
    xTo < exTo ?
    pdf.xfx(to,scale,xTo) :
    ((1.-xTo)/(1.-exTo)) * pdf.xfx(to,scale,exTo);

  if ( toPDF < 1e-8 )
    toPDF = 0.0;

  if ( toPDF == 0.0 || fromPDF == 0.0 )
    return 0.0;

  return toPDF / fromPDF;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PDFRatio::persistentOutput(PersistentOStream & os) const {
  os << theValenceExtrapolation << theSeaExtrapolation
     << ounit(theFreezingScale,GeV);
}

void PDFRatio::persistentInput(PersistentIStream & is, int) {
  is >> theValenceExtrapolation >> theSeaExtrapolation
     >> iunit(theFreezingScale,GeV);
}

ClassDescription<PDFRatio> PDFRatio::initPDFRatio;
// Definition of the static class description member.

void PDFRatio::Init() {

  static ClassDocumentation<PDFRatio> documentation
    ("PDFRatio implements numerically stable pdf ratios.");


  static Parameter<PDFRatio,double> interfaceValenceExtrapolation
    ("ValenceExtrapolation",
     "The x from which on extrapolation should be done for valence partons.",
     &PDFRatio::theValenceExtrapolation, 0.7, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<PDFRatio,double> interfaceSeaExtrapolation
    ("SeaExtrapolation",
     "The x from which on extrapolation should be done for valence partons.",
     &PDFRatio::theSeaExtrapolation, 0.6, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<PDFRatio,Energy> interfaceFreezingScale
    ("FreezingScale",
     "The scale below which the PDFs are frozen.",
     &PDFRatio::theFreezingScale, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

