// -*- C++ -*-
//
// MatchboxInsertionOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

using namespace Herwig;

MatchboxInsertionOperator::MatchboxInsertionOperator() 
  : HandlerBase(),
    theUseDR(false), theUseCS(false), theUseBDK(false), theUseExpanded(false) {}

MatchboxInsertionOperator::~MatchboxInsertionOperator() {}

void MatchboxInsertionOperator::setBorn(Ptr<MatchboxMEBase>::tptr me) { theLastBorn = me; }

Ptr<MatchboxMEBase>::tptr MatchboxInsertionOperator::lastBorn() { return theLastBorn; }

Ptr<MatchboxMEBase>::tcptr MatchboxInsertionOperator::lastBorn() const { return theLastBorn; }

CrossSection MatchboxInsertionOperator::dSigHatDR() const {
  return
    sqr(hbarc) * me2() *
    lastBorn()->lastXComb().jacobian() * 
    lastMEPDFWeight() /
    (2.*lastSHat());
}

void MatchboxInsertionOperator::additionalKinematics(const double * r) {
  if ( nDimAdditional() ) {
    additionalRandomNumbers.resize(nDimAdditional());
    copy(r,r+nDimAdditional(),additionalRandomNumbers.begin());
  }
}

void MatchboxInsertionOperator::rebind(const TranslationMap & trans) {
  theLastBorn = trans.translate(theLastBorn);
  HandlerBase::rebind(trans);
}

IVector MatchboxInsertionOperator::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.push_back(theLastBorn);
  return ret;
}

void MatchboxInsertionOperator::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theLastBorn
     << theUseDR << theUseCS << theUseBDK << theUseExpanded;
}

void MatchboxInsertionOperator::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theLastBorn
     >> theUseDR >> theUseCS >> theUseBDK >> theUseExpanded;
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
describeMatchboxInsertionOperator("Herwig::MatchboxInsertionOperator", "HwMatchbox.so");

