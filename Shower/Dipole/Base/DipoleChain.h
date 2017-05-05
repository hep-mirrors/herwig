// -*- C++ -*-
//
// DipoleChain.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleChain_H
#define HERWIG_DipoleChain_H
//
// This is the declaration of the DipoleChain class.
//

#include "Dipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief The DipoleChain class is used by the dipole shower to
 * represent a chain of dipoles.
 *
 */
class DipoleChain {

public:

  /**
   * Default constructor
   */
  DipoleChain();

  /**
   * Return true, if this chain is circular.
   */
  bool circular () const;

  /*
   * Return true, if the dipole referred to
   * has a left neighbour
   */
  bool hasLeftNeighbour(list<Dipole>::const_iterator dc) const;

  /*
   * Return a reference to the left neighbour,
   * if existing
   */
  Dipole& leftNeighbour(list<Dipole>::iterator dc);

  /*
   * Return a const reference to the left neighbour,
   * if existing
   */
  const Dipole& leftNeighbour(list<Dipole>::const_iterator dc) const;

  /**
   * Return an iterator to the left neighbour
   */
  list<Dipole>::iterator leftNeighbourIterator(list<Dipole>::iterator dc);

  /*
   * Return true, if the dipole referred to
   * has a right neighbour
   */
  bool hasRightNeighbour (list<Dipole>::const_iterator dc) const;

  /*
   * Return a reference to the right neighbour,
   * if existing
   */
  Dipole& rightNeighbour (list<Dipole>::iterator dc);

  /*
   * Return a const reference to the right neighbour,
   * if existing
   */
  const Dipole& rightNeighbour (list<Dipole>::const_iterator dc) const;

  /**
   * Return an iterator to the right neighbour
   */
  list<Dipole>::iterator rightNeighbourIterator(list<Dipole>::iterator dc);

public:

  /**
   * Access the dipole list
   */
  list<Dipole>& dipoles() { return theDipoles; }

  /**
   * Return the dipole list
   */
  const list<Dipole>& dipoles() const { return theDipoles; }

  /**
   * Check for gg single dipole
   */
  void check();

public:

  /*
   * Insert the given splitting; if this contains a chain-breakup emission and
   * the chain is circular, reshuffle the chain to make it non-circular; if it is
   * already non-circular return the iterator starting the new chain. If no
   * splitting is needed return the end iterator of the dipole list.
   * Set the iterators pointing to the children dipoles.
   */
  list<Dipole>::iterator insertSplitting(list<Dipole>::iterator emittingDipole,
					 pair<Dipole,Dipole> children,
					 pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators);

  /**
   * Indicate a change in the given dipole.
   */
  void updateDipole(list<Dipole>::iterator dip);

public:

  /**
   * Put information to ostream
   */
  void print(ostream&) const;

private:

  /**
   * The dipoles contained in this chain
   */
  list<Dipole> theDipoles;

  /**
   * Switch on special treatment for
   * gg single dipole
   */
  bool ggSingleDipole;

};

inline ostream& operator << (ostream& os, const DipoleChain& di) {
  di.print(os);
  return os;
}

}

#endif /* HERWIG_DipoleChain_H */
