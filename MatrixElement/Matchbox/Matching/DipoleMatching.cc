// -*- C++ -*-
//
// DipoleMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMatching class.
//

#include "DipoleMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

DipoleMatching::DipoleMatching() {}

DipoleMatching::~DipoleMatching() {}

IBPtr DipoleMatching::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMatching::fullclone() const {
  return new_ptr(*this);
}

CrossSection DipoleMatching::dSigHatDR() const {
  double pdfFactor = 1.;
  if ( showerScalesInSubtraction() ) {
    double bornPDF = bornPDFWeight(showerScalesInSubtraction());
    double bornPDFHard = bornPDFWeight(false);
    pdfFactor = bornPDFHard / bornPDF;
  }
  return
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    realPDFWeight(showerScalesInSubtraction()) * pdfFactor *
    couplingWeight(showerScalesInSubtraction()) *
    dipole()->me2() /
    (2. * realXComb()->lastSHat());
}

double DipoleMatching::me2() const {
  throw Exception() << "Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::abortnow;
  return 0.;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleMatching::persistentOutput(PersistentOStream &) const {}

void DipoleMatching::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMatching,Herwig::ShowerApproximation>
  describeHerwigDipoleMatching("Herwig::DipoleMatching", "HwMatchbox.so");

void DipoleMatching::Init() {

  static ClassDocumentation<DipoleMatching> documentation
    ("DipoleMatching implements naive NLO matching with the dipole shower.");

}

