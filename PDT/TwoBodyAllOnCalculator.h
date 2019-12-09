// -*- C++ -*-
//
// TwoBodyAllOnCalculator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoBodyAllOnCalculator_H
#define HERWIG_TwoBodyAllOnCalculator_H
// This is the declaration of the TwoBodyAllOnCalculator class.

#include "WidthCalculatorBase.h"
#include "TwoBodyAllOnCalculator.fh"
#include "GenericWidthGenerator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>TwoBodyAllOnCalculator</code> class is a wrapped around the the 
 *  simple two body decay matrix elements in the <code>GenericWidthGenerator</code>
 *  class and is designed to allow these matrix elements to be integrated
 *  if the external particles can be off-shell.
 *
 * @see TwoBodyAllOnCalculator
 * 
 */
class TwoBodyAllOnCalculator: public WidthCalculatorBase {

public:

  /**
   * The GenericWidthGenerator class is a friend to allow easier access for the
   * integration of the two body partial widths.
   */
  friend class GenericWidthGenerator;

public:

  /**
   * Constructor.
   * @param inwidth Pointer to the  GenericWidthGenerator class.
   * @param imode The mode in the GenericWidthGenerator class we are integrating
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   */
  TwoBodyAllOnCalculator(tGenericWidthGeneratorPtr inwidth,int imode,
			 Energy m1,Energy m2)
    : _mode(imode),_mass1(m1),_mass2(m2),_widthgen(inwidth)
  {}

  /**
   * member to calculate the partial width.
   * @param scale The mass squared for the decaying particle.
   * @return The partial width.
   */
  Energy partialWidth(Energy2 scale) const;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @param mass The new value.
   * @return The mass required.
   */
  void resetMass(int imass,Energy mass) {
    if(imass==1)      _mass1=mass;
    else if(imass==2) _mass2=mass;
    else throw Exception() << "Unknown particle in " 
			   << "TwoBodyAllOnCalculator::resetMass()"
			   << Exception::runerror;
  }

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  Energy getMass(const int imass) const {
    if(imass==1)      return _mass1;
    else if(imass==2) return _mass2;
    else throw Exception() << "Unknown particle in " 
			   << "TwoBodyAllOnCalculator::getMass()"
			   << Exception::runerror;
  }

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  Energy otherMass(const int imass) const {
    if(imass==1)      return _mass2;
    else if(imass==2) return _mass1;
    else throw Exception() << "Unknown particle in " 
			   << "TwoBodyAllOnCalculator::otherMass()"
			   << Exception::runerror;
  }

private:

  /**
   * Private and non-existent assignment operator.
   */
  TwoBodyAllOnCalculator & operator=(const TwoBodyAllOnCalculator &) = delete;

private:

  /**
   * the mode
   */
  int _mode;

  /**
   * Mass of the first particle.
   */
  Energy _mass1;

  /**
   * Mass of the second particle.
   */
  Energy _mass2;

  /**
   * the width generator
   */
  GenericWidthGeneratorPtr _widthgen;

};

}

#endif /* HERWIG_TwoBodyAllOnCalculator_H */
