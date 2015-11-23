// -*- C++ -*-
//
// MatchboxInsertionOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxInsertionOperator class.
//

#include "MatchboxInsertionOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

Ptr<MatchboxFactory>::tptr MatchboxInsertionOperator::factory() const { return theFactory; }

void MatchboxInsertionOperator::factory(Ptr<MatchboxFactory>::tptr f) { theFactory = f; }

MatchboxInsertionOperator::MatchboxInsertionOperator() 
  : HandlerBase(),
    theUseDRbar(false),
    theUseDR(false), theUseCS(false), theUseBDK(false), theUseExpanded(false) {}

MatchboxInsertionOperator::~MatchboxInsertionOperator() {}

void MatchboxInsertionOperator::cloneDependencies(const std::string&) {}

CrossSection MatchboxInsertionOperator::dSigHatDR() const {
  return
    sqr(hbarc) * me2() *
    lastBorn()->lastXComb().jacobian() * 
    lastMEPDFWeight() /
    (2.*lastSHat());
}

void MatchboxInsertionOperator::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theFactory << theUseDRbar
     << theUseDR << theUseCS << theUseBDK << theUseExpanded;
}

void MatchboxInsertionOperator::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theFactory >> theUseDRbar
     >> theUseDR >> theUseCS >> theUseBDK >> theUseExpanded;
  lastMatchboxXComb(theLastXComb);
}


void MatchboxInsertionOperator::Init() {

  static ClassDocumentation<MatchboxInsertionOperator> documentation
    ("MatchboxInsertionOperator is the base class for insertion operators");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxInsertionOperator,HandlerBase>
describeMatchboxInsertionOperator("Herwig::MatchboxInsertionOperator", "Herwig.so");

