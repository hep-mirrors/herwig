// -*- C++ -*-
//
// DipoleChain.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleChain class.
//

#include "DipoleChain.h"
#include "Herwig/DipoleShower/Utility/DipolePartonSplitter.h"

#include <boost/utility.hpp>

using namespace Herwig;

DipoleChain::DipoleChain() 
  : ggSingleDipole(false) {}

bool DipoleChain::circular () const { 
  return 
    (theDipoles.front().leftParticle() == 
     theDipoles.back().rightParticle());
}

bool DipoleChain::hasLeftNeighbour (list<Dipole>::const_iterator dc) const {
  if ( dc == dipoles().begin() )
    return circular();
  return true;
}

Dipole& DipoleChain::leftNeighbour (list<Dipole>::iterator dc) {
  assert(hasLeftNeighbour(dc));
  if ( dc == dipoles().begin() )
    return dipoles().back();
  return *(--dc);
}

const Dipole& DipoleChain::leftNeighbour (list<Dipole>::const_iterator dc) const {
  assert(hasLeftNeighbour(dc));
  if ( dc == dipoles().begin() )
    return dipoles().back();
  return *(--dc);
}

list<Dipole>::iterator DipoleChain::leftNeighbourIterator(list<Dipole>::iterator dc) {
  assert(hasLeftNeighbour(dc));
  if ( dc == dipoles().begin() )
    return --dipoles().end();
  return --dc;
}

bool DipoleChain::hasRightNeighbour (list<Dipole>::const_iterator dc) const {
  if (dc == --dipoles().end())
    return circular();
  return true;
}

Dipole& DipoleChain::rightNeighbour (list<Dipole>::iterator dc) {
  assert(hasRightNeighbour(dc));
  if ( dc == --dipoles().end() )
    return dipoles().front();
  return *(++dc);
}

const Dipole& DipoleChain::rightNeighbour (list<Dipole>::const_iterator dc) const {
  assert(hasRightNeighbour(dc));
  if ( dc == --dipoles().end() )
    return dipoles().front();
  return *(++dc);
}

list<Dipole>::iterator DipoleChain::rightNeighbourIterator(list<Dipole>::iterator dc) {
  assert(hasRightNeighbour(dc));
  if ( dc == --dipoles().end() )
    return dipoles().begin();
  return ++dc;
}

void DipoleChain::check() {
  if ( theDipoles.begin() == boost::prior(theDipoles.end()) ) {
    if ( theDipoles.front().leftParticle()->hasColour() &&
	 theDipoles.front().leftParticle()->hasAntiColour() ) {
      assert(theDipoles.front().rightParticle()->hasColour() &&
	     theDipoles.front().rightParticle()->hasAntiColour());
      ggSingleDipole = true;
    }
  }
}

list<Dipole>::iterator DipoleChain::insertSplitting(list<Dipole>::iterator emittingDipole,
						    pair<Dipole,Dipole> children,
						    pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators) {

  assert(DipolePartonSplitter::colourConnected(children.first.leftParticle(),children.first.rightParticle()) ||
	 DipolePartonSplitter::colourConnected(children.second.leftParticle(),children.second.rightParticle()));

  bool was_circular = circular();

  if (hasLeftNeighbour(emittingDipole)) {

    list<Dipole>::iterator theLeftNeighbour =
      leftNeighbourIterator(emittingDipole);

    theLeftNeighbour->rightParticle(children.first.leftParticle());
    if ( children.first.leftParticle()->scale() < sqr(theLeftNeighbour->rightScale()) )
      theLeftNeighbour->rightScale(sqrt(children.first.leftParticle()->scale()));
    theLeftNeighbour->rightPDF(children.first.leftPDF());
    theLeftNeighbour->rightFraction(children.first.leftFraction());

    theLeftNeighbour->update();

  }

  if (hasRightNeighbour(emittingDipole)) {

    list<Dipole>::iterator theRightNeighbour =
      rightNeighbourIterator(emittingDipole);

    theRightNeighbour->leftParticle(children.second.rightParticle());
    if ( children.second.rightParticle()->scale() < sqr(theRightNeighbour->leftScale()) )
      theRightNeighbour->leftScale(sqrt(children.second.rightParticle()->scale()));
    theRightNeighbour->leftPDF(children.second.rightPDF());
    theRightNeighbour->leftFraction(children.second.rightFraction());

    theRightNeighbour->update();

  }

  if (DipolePartonSplitter::colourConnected(children.first.leftParticle(),children.first.rightParticle()) &&
      DipolePartonSplitter::colourConnected(children.second.leftParticle(),children.second.rightParticle())) {

    // nothing special to do, just replace the emitting dipole
    // by the right one and insert the left one before it

    *emittingDipole = children.second;

    childIterators.second = emittingDipole;
    childIterators.first = dipoles().insert(emittingDipole,children.first);

    if ( ggSingleDipole ) {
      ggSingleDipole = false;
      Dipole miss;
      miss.leftParticle(dipoles().back().rightParticle());
      miss.rightParticle(dipoles().front().leftParticle());
      miss.leftScale(dipoles().back().rightScale());
      miss.rightScale(dipoles().front().leftScale());
      miss.leftPDF(dipoles().back().rightPDF());
      miss.rightPDF(dipoles().front().leftPDF());
      miss.leftFraction(dipoles().back().rightFraction());
      miss.rightFraction(dipoles().front().leftFraction());
      miss.update();
      dipoles().push_back(miss);
    }

    return dipoles().end();

  }

  if (!DipolePartonSplitter::colourConnected(children.first.leftParticle(),children.first.rightParticle())) {

    if ( !was_circular && !ggSingleDipole ) {
      *emittingDipole = children.second;
      childIterators.second = emittingDipole;
      assert(emittingDipole != dipoles().begin());
      childIterators.first = boost::prior(emittingDipole);
      return emittingDipole;
    }

    *emittingDipole = children.second;

    if ( ggSingleDipole ) {
      ggSingleDipole = false;
      Dipole miss;
      miss.leftParticle(children.second.rightParticle());
      miss.rightParticle(children.first.leftParticle());
      miss.leftScale(children.second.rightScale());
      miss.rightScale(children.first.leftScale());
      miss.leftPDF(children.second.rightPDF());
      miss.rightPDF(children.first.leftPDF());
      miss.leftFraction(children.second.rightFraction());
      miss.rightFraction(children.first.leftFraction());
      miss.update();
      dipoles().push_back(miss);
      childIterators.first = dipoles().begin();
      childIterators.second = boost::prior(dipoles().end());
      return dipoles().end();
    }

    childIterators.second = emittingDipole;
    if ( emittingDipole == dipoles().begin() )
      childIterators.first = --dipoles().end();
    else
      childIterators.first = boost::prior(emittingDipole);

    if ( emittingDipole == dipoles().begin() )
      return dipoles().end();

    dipoles().splice(dipoles().begin(),dipoles(),emittingDipole,dipoles().end());

    // explicitly fix iterators in case the splice implementation
    // at hand does invalidate iterators (the SGI docu says, it doesn't,
    // but it seems that this behaviour is not part of the standard)
    childIterators.second = dipoles().begin();
    childIterators.first = --dipoles().end();
    
    return dipoles().end();

  }

  if (!DipolePartonSplitter::colourConnected(children.second.leftParticle(),children.second.rightParticle())) {

    if ( !was_circular && !ggSingleDipole ) {
      *emittingDipole = children.first;
      childIterators.first = emittingDipole;
      assert(emittingDipole != --dipoles().end());
      childIterators.second = boost::next(emittingDipole);
      return boost::next(emittingDipole);
    }

    *emittingDipole = children.first;

    if ( ggSingleDipole ) {
      ggSingleDipole = false;
      Dipole miss;
      miss.leftParticle(children.second.rightParticle());
      miss.rightParticle(children.first.leftParticle());
      miss.leftScale(children.second.rightScale());
      miss.rightScale(children.first.leftScale());
      miss.leftPDF(children.second.rightPDF());
      miss.rightPDF(children.first.leftPDF());
      miss.leftFraction(children.second.rightFraction());
      miss.rightFraction(children.first.leftFraction());
      miss.update();
      dipoles().push_front(miss);
      childIterators.first = dipoles().begin();
      childIterators.second = boost::prior(dipoles().end());
      return dipoles().end();
    }

    childIterators.first = emittingDipole;
    if ( emittingDipole == --dipoles().end() )
      childIterators.second = dipoles().begin();
    else
      childIterators.second = boost::next(emittingDipole);

    if ( emittingDipole == --dipoles().end() )
      return dipoles().end();

    dipoles().splice(dipoles().begin(),dipoles(),boost::next(emittingDipole),dipoles().end());

    // explicitly fix iterators in case the splice implementation
    // at hand does invalidate iterators (the SGI docu says, it doesn't,
    // but it seems that this behaviour is not part of the standard)
    childIterators.first = dipoles().begin();
    childIterators.second = --dipoles().end();
    
    return dipoles().end();

  }

  return dipoles().end();  

}

void DipoleChain::updateDipole(list<Dipole>::iterator dip) {

  dip->update();

  if (hasLeftNeighbour(dip)) {

    list<Dipole>::iterator theLeftNeighbour =
      leftNeighbourIterator(dip);

    theLeftNeighbour->rightParticle(dip->leftParticle());
    theLeftNeighbour->rightPDF(dip->leftPDF());
    theLeftNeighbour->rightFraction(dip->leftFraction());

    theLeftNeighbour->update();

  }

  if (hasRightNeighbour(dip)) {

    list<Dipole>::iterator theRightNeighbour =
      rightNeighbourIterator(dip);

    theRightNeighbour->leftParticle(dip->rightParticle());
    theRightNeighbour->leftPDF(dip->rightPDF());
    theRightNeighbour->leftFraction(dip->rightFraction());

    theRightNeighbour->update();

  }

}

void DipoleChain::print(ostream& os) const {

  os << "--- DipoleChain ----------------------------------------------------------------\n";

  if ( theDipoles.empty() ) {
    os << "  ***  This DipoleChain is empty.  ***\n";
  } else {

    os << " " << (!circular() ? "non-" : "") << "circular with "
       << theDipoles.size() << " dipoles\n";

    for (list<Dipole>::const_iterator dit = theDipoles.begin();
	 dit != theDipoles.end(); ++dit) {
      os << (*dit);
    }

  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}
