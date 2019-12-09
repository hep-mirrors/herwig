// -*- C++ -*-
//
// MatchboxTopMTScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxTopMTScale class.
//

#include "MatchboxTopMTScale.h"
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

MatchboxTopMTScale::MatchboxTopMTScale() :
  theShowerScaleMode(1) {}

MatchboxTopMTScale::~MatchboxTopMTScale() {}

IBPtr MatchboxTopMTScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxTopMTScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxTopMTScale::renormalizationScale() const {
  
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
    throw Exception() << "MatchboxTopMTScale: Could not find a top-antitop-pair in the final state!\n"
		      << Exception::runerror;
  }
  // cerr << " sqrt(TopMTScale)  = "
  //      << sqrt(meMomenta()[top].mt2()+meMomenta()[antitop].mt2())/GeV
  //      << "\n" << flush;
  return 0.5*(meMomenta()[top].mt2()+meMomenta()[antitop].mt2());

}

Energy2 MatchboxTopMTScale::factorizationScale() const {
  return(renormalizationScale());
}


Energy2 MatchboxTopMTScale::showerScale() const {
  
  if ( theShowerScaleMode == 1 )
    return factorizationScale();
  
  else {
    assert(theShowerScaleMode == 2);
    size_t k = 2;
    int top = -1;
    int antitop = -1;
    int emission = -1;
  
    if ( mePartonData().size() > 5 )
      assert(false && "MatchboxTopMTScale.cc: mePartonData().size() > 5.\n");

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
	  assert(false && "More than one emission in MatchboxTopMTScale.\n");
      }
    
      k++;
    }

    if ( top < 2 || antitop < 2 ){
      throw Exception() << "MatchboxTopMTScale: Could not find a top-antitop-pair in the final state!\n"
			<< Exception::runerror;
    }
  

    if ( mePartonData().size() > 4 && emission < 2 ){
      throw Exception() << "MatchboxTopMTScale: Could not find an emission in a state where mePartonData().size() > 4.\n"
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


void MatchboxTopMTScale::persistentOutput(PersistentOStream & os) const {
  os << theShowerScaleMode;
}

void MatchboxTopMTScale::persistentInput(PersistentIStream & is, int) {
  is >> theShowerScaleMode;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxTopMTScale,MatchboxScaleChoice>
  describeHerwigMatchboxTopMTScale("Herwig::MatchboxTopMTScale", "HwMatchboxScales.so");

void MatchboxTopMTScale::Init() {

  static ClassDocumentation<MatchboxTopMTScale> documentation
    ("MatchboxTopMTScale implements the quadratic sum of the transverse masses of the top and antitop quark as a scale choice.");


  static Switch<MatchboxTopMTScale, unsigned int> interfaceShowerScaleMode
    ("ShowerScaleMode",
     "Choose the definition of the shower hard scale.",
     &MatchboxTopMTScale::theShowerScaleMode, 1, false, false);
  static SwitchOption FactorizationScale
    (interfaceShowerScaleMode,"FactorizationScale","Use the factorization scale.", 1);
  static SwitchOption MeanMT2
    (interfaceShowerScaleMode,"MeanMT2","Use the mean squared transverse mass of the outgoing particles.", 2);

}

