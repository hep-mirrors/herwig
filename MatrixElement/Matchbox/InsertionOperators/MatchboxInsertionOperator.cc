// -*- C++ -*-
//
// MatchboxInsertionOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

Ptr<MatchboxFactory>::tptr MatchboxInsertionOperator::factory() const {
  return MatchboxFactory::currentFactory();
}

MatchboxInsertionOperator::MatchboxInsertionOperator() 
  : HandlerBase() {}

MatchboxInsertionOperator::~MatchboxInsertionOperator() {}

void MatchboxInsertionOperator::cloneDependencies(const std::string&) {}

CrossSection MatchboxInsertionOperator::dSigHatDR() const {
  return
    sqr(hbarc) * me2() *
    lastBorn()->lastXComb().jacobian() * 
    lastMEPDFWeight() /
    (2.*lastSHat());
}

Ptr<MatchboxMEBase>::tptr MatchboxInsertionOperator::lastBorn() const {
    if(lastMatchboxXComb()->matchboxME()) return lastMatchboxXComb()->matchboxME();
     else{
       assert(lastMatchboxXComb()->subtractionDipole());
       return lastMatchboxXComb()->subtractionDipole()->underlyingBornME();
     }
}



bool MatchboxInsertionOperator::isDRbar() const { return lastBorn()->isDRbar(); }

bool MatchboxInsertionOperator::isDR() const { return lastBorn()->isDR(); }

bool MatchboxInsertionOperator::isCS() const { return lastBorn()->isCS(); }

bool MatchboxInsertionOperator::isBDK() const { return lastBorn()->isBDK(); }

bool MatchboxInsertionOperator::isExpanded() const { return lastBorn()->isExpanded(); }






void MatchboxInsertionOperator::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb;
}

void MatchboxInsertionOperator::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb ;
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

