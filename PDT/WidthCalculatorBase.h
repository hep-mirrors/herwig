// -*- C++ -*-
//
// WidthCalculatorBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WidthCalculatorBase_H
#define HERWIG_WidthCalculatorBase_H
//
// This is the declaration of the WidthCalculatorBase class.
//
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Pointer/ReferenceCounted.h"
#include "WidthCalculatorBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *  The WidthCalculatorBase class is a base class to be used
 *  by classes which calculate partial widths for the running width.
 *
 * @see DecayIntegrator
 * @see GenericWidthGenerator
 * 
 */
class WidthCalculatorBase: public Pointer::ReferenceCounted {

public:

  /**
   *  Destructor
   */
  virtual ~WidthCalculatorBase();

  /**
   * Calculate the partial width. This must be implemented in classes inheriting from
   * this one.
   * @param scale The mass squared of the decaying particle.
   * @return The partial width.
   */
  virtual Energy partialWidth(Energy2 scale) const =0;

  /**
   * Reset the mass of a particle (used to integrate over the mass.) This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass to be reset.
   * @param mass The new mass.
   */
  virtual void resetMass(int imass,Energy mass) =0;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  virtual Energy getMass(const int imass) const= 0;

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  virtual Energy otherMass(const int imass) const=0;

private:

  /**
   * Private and non-existent assignment operator.
   */
  WidthCalculatorBase & operator=(const WidthCalculatorBase &) = delete;

};
}

#endif /* HERWIG_WidthCalculatorBase_H */
