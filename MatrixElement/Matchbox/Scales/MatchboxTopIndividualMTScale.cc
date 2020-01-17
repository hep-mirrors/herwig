// -*- C++ -*-
//
// MatchboxTopIndividualMTScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxTopIndividualMTScale class.
//

#include "MatchboxTopIndividualMTScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxTopIndividualMTScale::MatchboxTopIndividualMTScale() :
  theShowerScaleMode(1), theFactor(1.), theTopOrAntitop(1) {}

MatchboxTopIndividualMTScale::~MatchboxTopIndividualMTScale() {}

IBPtr MatchboxTopIndividualMTScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxTopIndividualMTScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxTopIndividualMTScale::renormalizationScale() const {
  
  size_t k = 2;
  int top = -1;
  int antitop = -1;
  
  while ( (top == -1 || antitop == -1) && k < mePartonData().size() ){
    if ( mePartonData()[k]->id() == 6 ) {
      if ( top < 0 )
	top = k;
      else
	assert(false);
    } else if ( mePartonData()[k]->id() == -6 ) {
      if ( antitop < 0 )
	antitop = k;
      else
	assert(false);
    }
    k++;
  }

  if ( top < 2 || antitop < 2 ){
    throw Exception() << "MatchboxTopIndividualMTScale: Could not find a top-antitop-pair in the final state!\n"
		      << Exception::runerror;
  }

  Energy scale = ZERO;
  
  // Not using .mt() as this is signed and not what we want.

  if ( theTopOrAntitop == 1 )
    scale = sqrt(meMomenta()[top].mt2());
  else {
    assert ( theTopOrAntitop == 2);
    scale = sqrt(meMomenta()[antitop].mt2());
  }
  
  return sqr(scale*theFactor);

}

Energy2 MatchboxTopIndividualMTScale::factorizationScale() const {
  return(renormalizationScale());
}

Energy2 MatchboxTopIndividualMTScale::showerScale() const {
  
  if ( theShowerScaleMode == 1 )
    return factorizationScale();
  
  else {
    assert(theShowerScaleMode == 2);
    size_t k = 2;
    int top = -1;
    int antitop = -1;
    int emission = -1;
  
    if ( mePartonData().size() > 5 )
      assert(false && "MatchboxTopIndividualMTScale.cc: mePartonData().size() > 5.\n");

    while ( (top == -1 || antitop == -1 || emission == -1) && k < mePartonData().size() ){
      if ( mePartonData()[k]->id() == 6 ) {
	if ( top < 0 )
	  top = k;
	else
	  assert(false);
      } else if ( mePartonData()[k]->id() == -6 ) {
	if ( antitop < 0 )
	  antitop = k;
	else
	  assert(false);
      }
      else if ( abs(mePartonData()[k]->id()) < 6 || mePartonData()[k]->id() == 21 ) {
	if ( emission < 0 )
	  emission = k;
	else
	  assert(false && "More than one emission in MatchboxTopIndividualMTScale.\n");
      }
    
      k++;
    }

    if ( top < 2 || antitop < 2 ){
      throw Exception() << "MatchboxTopIndividualMTScale: Could not find a top-antitop-pair in the final state!\n"
			<< Exception::runerror;
    }
  

    if ( mePartonData().size() > 4 && emission < 2 ){
      throw Exception() << "MatchboxTopIndividualMTScale: Could not find an emission in a state where mePartonData().size() > 4.\n"
			<< Exception::runerror;
    }

    if ( emission < 2 )
      return (1./2.)*(meMomenta()[top].mt2()+meMomenta()[antitop].mt2());
  
    else
      return (1./3.)*(meMomenta()[top].mt2()+meMomenta()[antitop].mt2() + meMomenta()[emission].mt2());
  }
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxTopIndividualMTScale::persistentOutput(PersistentOStream & os) const {
  os << theShowerScaleMode << theFactor << theTopOrAntitop;
}

void MatchboxTopIndividualMTScale::persistentInput(PersistentIStream & is, int) {
  is >> theShowerScaleMode >> theFactor >> theTopOrAntitop;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxTopIndividualMTScale,MatchboxScaleChoice>
  describeHerwigMatchboxTopIndividualMTScale("Herwig::MatchboxTopIndividualMTScale", "HwMatchboxScales.so");

void MatchboxTopIndividualMTScale::Init() {

  static ClassDocumentation<MatchboxTopIndividualMTScale> documentation
    ("MatchboxTopIndividualMTScale implements the transverse mass of the top quark as the scale choice.");

  static Switch<MatchboxTopIndividualMTScale, unsigned int> interfaceShowerScaleMode
    ("ShowerScaleMode",
     "Choose the definition of the shower hard scale.",
     &MatchboxTopIndividualMTScale::theShowerScaleMode, 1, false, false);
  static SwitchOption FactorizationScale
    (interfaceShowerScaleMode,"FactorizationScale","Use the factorization scale.", 1);
  static SwitchOption MeanMT2
    (interfaceShowerScaleMode,"MeanMT2","Use the mean squared transverse mass of the outgoing particles.", 2);

  
  static Parameter<MatchboxTopIndividualMTScale,double> interfaceMultiplicationFactor
    ("MultiplicationFactor",
     "Set a multiplicative factor to include in the scale choice definition.",
     &MatchboxTopIndividualMTScale::theFactor, 1., 0., 0,
     false, false, Interface::lowerlim);


  static Switch<MatchboxTopIndividualMTScale, unsigned int> interfaceTopOrAntitop
    ("TopOrAntitop",
     "Choose which quark transverse mass to use.",
     &MatchboxTopIndividualMTScale::theTopOrAntitop, 1, false, false);
  static SwitchOption Top
    (interfaceTopOrAntitop,"Top","Use the top quark transverse mass.", 1);
  static SwitchOption Antitop
    (interfaceTopOrAntitop,"Antitop","Use the top antiquark transverse mass.", 2);
  
}

