// -*- C++ -*-
//
// EnhanceNLOContributions.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EnhanceNLOContributions class.
//

#include "EnhanceNLOContributions.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/SubtractedME.h"

using namespace Herwig;

EnhanceNLOContributions::EnhanceNLOContributions() 
  : theEnhanceInitialPointsVirtual(1.0), theOversamplingFactorVirtual(1.0),
    theEnhanceInitialPointsReal(1.0), theOversamplingFactorReal(1.0) {}

EnhanceNLOContributions::~EnhanceNLOContributions() {}

IBPtr EnhanceNLOContributions::clone() const {
  return new_ptr(*this);
}

IBPtr EnhanceNLOContributions::fullclone() const {
  return new_ptr(*this);
}

double EnhanceNLOContributions::enhanceInitialPoints(const StandardXComb& xc) const {
  Ptr<SubtractedME>::tptr subme = dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(xc.matrixElement());
  if ( subme ) {
    return theEnhanceInitialPointsReal;
  }
  Ptr<MatchboxMEBase>::tptr me = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(xc.matrixElement());
  if ( me )
    if ( me->oneLoopNoBorn() ) {
      return theEnhanceInitialPointsVirtual;
    }
  return 1.0;
}

double EnhanceNLOContributions::oversamplingFactor(const StandardXComb& xc) const {
  Ptr<SubtractedME>::tptr subme = dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(xc.matrixElement());
  if ( subme ) {
    return theOversamplingFactorReal;
  }
  Ptr<MatchboxMEBase>::tptr me = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(xc.matrixElement());
  if ( me )
    if ( me->oneLoopNoBorn() ) {
      return theOversamplingFactorVirtual;
    }
  return 1.0;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void EnhanceNLOContributions::persistentOutput(PersistentOStream & os) const {
  os << theEnhanceInitialPointsVirtual << theOversamplingFactorVirtual
     << theEnhanceInitialPointsReal << theOversamplingFactorReal;
}

void EnhanceNLOContributions::persistentInput(PersistentIStream & is, int) {
  is >> theEnhanceInitialPointsVirtual >> theOversamplingFactorVirtual
     >> theEnhanceInitialPointsReal >> theOversamplingFactorReal;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<EnhanceNLOContributions,Herwig::SamplingBias>
  describeHerwigEnhanceNLOContributions("Herwig::EnhanceNLOContributions", 
				       "HwExsample2.so HwExsampleMatchbox.so");

void EnhanceNLOContributions::Init() {

  static ClassDocumentation<EnhanceNLOContributions> documentation
    ("Enhance NLO virtual or real emission contributions.");


  static Parameter<EnhanceNLOContributions,double> interfaceEnhanceInitialPointsVirtual
    ("EnhanceInitialPointsVirtual",
     "Enhance or supress the number initial points to integrate virtual contributions.",
     &EnhanceNLOContributions::theEnhanceInitialPointsVirtual, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<EnhanceNLOContributions,double> interfaceOversamplingFactorVirtual
    ("OversamplingFactorVirtual",
     "Enhance or supress the sampling of virtual contributions.",
     &EnhanceNLOContributions::theOversamplingFactorVirtual, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<EnhanceNLOContributions,double> interfaceEnhanceInitialPointsReal
    ("EnhanceInitialPointsReal",
     "Enhance or supress the number initial points to integrate real contributions.",
     &EnhanceNLOContributions::theEnhanceInitialPointsReal, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<EnhanceNLOContributions,double> interfaceOversamplingFactorReal
    ("OversamplingFactorReal",
     "Enhance or supress the sampling of real contributions.",
     &EnhanceNLOContributions::theOversamplingFactorReal, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

}

