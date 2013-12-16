// -*- C++ -*-
//
// NLOJetRandomWrapper.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_NLOJetRandomWrapper_H
#define HERWIG_NLOJetRandomWrapper_H

#include <cassert>
#include "nlo++/bits/hep-rng.h"

namespace Herwig {

/**
 * \ingroup Matchbox 
 * \author Simon Platzer 
 * \brief Wrap a nlo::random_generator around a vector of doubles
 */
class NLOJetRandomWrapper : public nlo::random_generator {

public:

  /**
   * Default constructor
   */
  NLOJetRandomWrapper()
    : random_generator(0,0), theNumbers(0) {}

  /**
   * Set the random numbers to be used in the next subsequent calls to
   * get_double()
   */
  void numbers(const double * newNumbers) { theNumbers = newNumbers; }

  /**
   * Set the seed -- irrelevant for this implementation.
   */
  virtual void set(unsigned long int) {}

  /**
   * Return a random integer. As up to now unsupported; is this needed
   * somewhere?
   */
  virtual unsigned long int get() const { assert(false && "implementation missing"); return 0; }

  /**
   * Return a random double.
   */
  virtual double get_double() const { double r = theNumbers[0]; ++theNumbers; return r; }

private:

  /**
   * The random numbers to be used in the next subsequent calls to
   * get_double()
   */
  mutable const double * theNumbers;

};

}

#endif // HERWIG_NLOJetRandomWrapper_H
