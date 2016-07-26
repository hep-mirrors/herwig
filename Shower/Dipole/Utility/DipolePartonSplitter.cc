// -*- C++ -*-
//
// DipolePartonSplitter.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipolePartonSplitter class.
//

#include "DipolePartonSplitter.h"

#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ColourLine.h"

using namespace Herwig;

void DipolePartonSplitter::split(tPPtr parent, tPPtr firstChild, tPPtr secondChild,
				 bool initialState) {

  firstChild->colourInfo(new_ptr(ColourBase()));
  secondChild->colourInfo(new_ptr(ColourBase()));

  if (!initialState) {
    parent->addChild(firstChild);
    parent->addChild(secondChild);
  } else {

    tPVector parents = parent->parents();

    for (tPVector::const_iterator q = parents.begin();
	 q != parents.end(); ++q) {

      (**q).addChild(firstChild);
      (**q).abandonChild(parent);

    }

    firstChild->addChild(parent);
    firstChild->addChild(secondChild);

  }

}

void DipolePartonSplitter::split(tPPtr parent, tPPtr firstChild, tPPtr secondChild,
				 tPPtr ref, bool initialState) {

  split(parent,firstChild,secondChild,initialState);

  // final state splittings

  if (!initialState) {

    // triplet radiating octet
    if (parent->hasColour() &&
	!parent->hasAntiColour()) {

      assert(secondChild->hasColour() &&
	     secondChild->hasAntiColour() &&
	     firstChild->hasColour() &&
	     !firstChild->hasAntiColour());

      secondChild->incomingColour(parent);
      firstChild->antiColourNeighbour(secondChild);

      return;

    }

    // anti-triplet radiating octet
    if (!parent->hasColour() &&
	parent->hasAntiColour()) {

      assert(secondChild->hasColour() &&
	     secondChild->hasAntiColour() &&
	     !firstChild->hasColour() &&
	     firstChild->hasAntiColour());

      secondChild->incomingAntiColour(parent);
      firstChild->colourNeighbour(secondChild);

      return;

    }

    // octet radiating octet
    if (parent->hasColour() &&
	parent->hasAntiColour() &&
	secondChild->hasColour() &&
	secondChild->hasAntiColour()) {

      assert(firstChild->hasColour() &&
	     firstChild->hasAntiColour());

      // check wether we radiated from the colour or anticolour line

      if (parent->colourLine() == ref->antiColourLine()) {
	firstChild->incomingAntiColour(parent);
	secondChild->incomingColour(parent);
	firstChild->antiColourNeighbour(secondChild);
      } else {
	firstChild->incomingColour(parent);
	secondChild->incomingAntiColour(parent);
	firstChild->colourNeighbour(secondChild);
      }

      return;

    }

    // octet splitting into triplet x anti-triplet
    if (parent->hasColour() &&
	parent->hasAntiColour() &&
	( (secondChild->hasAntiColour() &&
	   !secondChild->hasColour()) || 
	  (!secondChild->hasAntiColour() &&
	   secondChild->hasColour()) )) {

      if (firstChild->hasColour()) {
	assert(!firstChild->hasAntiColour() &&
	       secondChild->hasAntiColour() &&
	       !secondChild->hasColour());

	firstChild->incomingColour(parent);
	secondChild->incomingAntiColour(parent);

	return;

      }

      if (firstChild->hasAntiColour()) {
	assert(!firstChild->hasColour() &&
	       secondChild->hasColour() &&
	       !secondChild->hasAntiColour());

	firstChild->incomingAntiColour(parent);
	secondChild->incomingColour(parent);

	return;

      }

    }

  }

  // initial state splittings

  if (initialState) {

    // triplet
    if (parent->hasColour() &&
	!parent->hasAntiColour()) {

      if (secondChild->hasColour() &&
	  secondChild->hasAntiColour()) {

	// radiating octet

	assert(firstChild->hasColour() &&
	       !firstChild->hasAntiColour());

	parent->antiColourNeighbour(secondChild);
	firstChild->outgoingColour(secondChild);

      } else {

	// 'radiating' anti-triplet

	assert(!secondChild->hasColour() &&
	       secondChild->hasAntiColour() &&
	       firstChild->hasColour() &&
	       firstChild->hasAntiColour());

	secondChild->incomingAntiColour(firstChild);
	parent->colourLine()->addColoured(firstChild);

      }

      return;

    }

    // anti-triplet
    if (!parent->hasColour() &&
	parent->hasAntiColour()) {

      if (secondChild->hasColour() &&
	  secondChild->hasAntiColour()) {

	// radiating octet

	assert(!firstChild->hasColour() &&
	       firstChild->hasAntiColour());

	parent->colourNeighbour(secondChild);
	firstChild->outgoingAntiColour(secondChild);

      } else {

	// 'radiating' triplet

	assert(secondChild->hasColour() &&
	       !secondChild->hasAntiColour() &&
	       firstChild->hasColour() &&
	       firstChild->hasAntiColour());

	secondChild->incomingColour(firstChild);
	parent->antiColourLine()->addAntiColoured(firstChild);

      }

      return;

    }

    // octet radiating octet

    if (parent->hasColour() &&
	parent->hasAntiColour() &&
	secondChild->hasColour() &&
	secondChild->hasAntiColour()) {

      assert(firstChild->hasColour() &&
	     firstChild->hasAntiColour());

      // check wether we radiated from the colour or anticolour line

      if (parent->colourLine() == ref->colourLine()) {
	parent->colourLine()->addAntiColoured(secondChild);
	parent->antiColourLine()->addAntiColoured(firstChild);
	secondChild->incomingColour(firstChild);
      } else {
	parent->antiColourLine()->addColoured(secondChild);
	parent->colourLine()->addColoured(firstChild);
	secondChild->incomingAntiColour(firstChild);
      }

      return;

    }

    // octet splitting into triplet x triplet
    if (parent->hasColour() &&
	parent->hasAntiColour() &&
	( (secondChild->hasAntiColour() &&
	   !secondChild->hasColour()) || 
	  (!secondChild->hasAntiColour() &&
	   secondChild->hasColour()) )) {

      if (firstChild->hasColour()) {
	assert(!firstChild->hasAntiColour() &&
	       secondChild->hasColour() &&
	       !secondChild->hasAntiColour());
	
	parent->colourLine()->addColoured(firstChild);
	parent->antiColourLine()->addColoured(secondChild);

	return;

      }

      if (firstChild->hasAntiColour()) {
	assert(!firstChild->hasColour() &&
	       secondChild->hasAntiColour() &&
	       !secondChild->hasColour());

	parent->antiColourLine()->addAntiColoured(firstChild);
	parent->colourLine()->addAntiColoured(secondChild);

	return;

      }

    }


  }

  // should never happen
  assert(false);

}

void DipolePartonSplitter::change(tPPtr parent, tPPtr child, bool initialState) {

  child->colourInfo(new_ptr(ColourBase()));

  if (parent->hasColour())
    parent->colourLine()->addColoured(child);

  if (parent->hasAntiColour())
    parent->antiColourLine()->addAntiColoured(child);

  if (!initialState) {
    parent->addChild(child);
  } else {


    tPVector parents = parent->parents();

    for (tPVector::const_iterator q = parents.begin();
	 q != parents.end(); ++q) {

      (**q).addChild(child);
      (**q).abandonChild(parent);

    }

    child->addChild(parent);

  }

}

bool DipolePartonSplitter::colourConnected(tcPPtr first, tcPPtr second) {
  
  if (!first->coloured() || !second->coloured())
    return false;

  if (first->colourInfo()->colourLine()) {

    if (second->colourInfo()->colourLine() &&
	first->colourInfo()->colourLine() == second->colourInfo()->colourLine())
      return true;

    if (second->colourInfo()->antiColourLine() &&
	first->colourInfo()->colourLine() == second->colourInfo()->antiColourLine())
      return true;

  }

  if (first->colourInfo()->antiColourLine()) {

    if (second->colourInfo()->colourLine() &&
	first->colourInfo()->antiColourLine() == second->colourInfo()->colourLine())
      return true;

    if (second->colourInfo()->antiColourLine() &&
	first->colourInfo()->antiColourLine() == second->colourInfo()->antiColourLine())
      return true;

  }

  return false;

}
								  
