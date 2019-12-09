// -*- C++ -*-
//
// ColourMatrixElementCorrection.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourMatrixElementCorrection class.
//

#include "ColourMatrixElementCorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"

using namespace Herwig;

ColourMatrixElementCorrection::ColourMatrixElementCorrection():lambda(1.0),negCMECScaling(1.0) {}

ColourMatrixElementCorrection::~ColourMatrixElementCorrection() {}

IBPtr ColourMatrixElementCorrection::clone() const {
  return new_ptr(*this);
}

IBPtr ColourMatrixElementCorrection::fullclone() const {
  return new_ptr(*this);
}

double ColourMatrixElementCorrection::cmec(const DipoleSplittingInfo& dsplit) const {
  // Do not calculate CMECs if the subleading Nc shower has stopped 
  if ( !(currentHandler()->continueSubleadingNc()) ) {
    return 1.;
  }
  const PPtr em = dsplit.emitter();
  const PPtr sp = dsplit.spectator();
  const tcPDPtr emis = dsplit.emissionData();
  // Get the dictionary from particle pointers to herwig particle numbers
  const map<PPtr,size_t>& theDictionary = currentHandler()->particleIndices();
  const cPDVector& theParticles = currentHandler()->particlesAfter();
  // Use the dictionary to find
  // - emitter index
  // - spectator index
  // - emission ParticleID
  const std::tuple<size_t,size_t,long> ikemission = std::make_tuple(theDictionary.at(em),
								    theDictionary.at(sp),
								    emis->id());
  double factor = currentHandler()->densityOperator().colourMatrixElementCorrection(ikemission,theParticles);

  return factor > 0.0 ? factor : negCMECScaling*factor;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ColourMatrixElementCorrection::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ColourMatrixElementCorrection::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ColourMatrixElementCorrection,Herwig::DipoleSplittingReweight>
  describeHerwigColourMatrixElementCorrection("Herwig::ColourMatrixElementCorrection", 
					      "HwDipoleShower.so");

void ColourMatrixElementCorrection::Init() {

  static ClassDocumentation<ColourMatrixElementCorrection> documentation
    ("There is no documentation for the ColourMatrixElementCorrection class");

}

