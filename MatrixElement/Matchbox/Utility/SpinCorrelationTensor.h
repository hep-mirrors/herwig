// -*- C++ -*-
//
// SpinCorrelationTensor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SpinCorrelationTensor_H
#define HERWIG_SpinCorrelationTensor_H

#include "ThePEG/Config/Unitsystem.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief SpinCorrelationTensor represents a spin correlation tensor
 * of the form
 *
 * c g^{\mu\nu} + p^\mu p^\nu / Q^2
 *
 */
class SpinCorrelationTensor {

  /**
   * The diagonal part
   */
  double theDiagonal;

  /**
   * The vector to define to off-diagonal part
   */
  Lorentz5Momentum theMomentum;

  /**
   * The scale for the off-diagonal part
   */
  Energy2 theScale;

public:

  /**
   * The standard constructor, giving diagonal
   * part, and vector for off-diagonal components
   * and the scale.
   */
  SpinCorrelationTensor(double diag,
			const Lorentz5Momentum& p,
			Energy2 sc) 
    : theDiagonal(diag), theMomentum(p), theScale(sc) {}

public:

  /**
   * Return the diagonal part
   */
  double diagonal() const { return theDiagonal; }

  /**
   * Return the off-diagonal vector part.
   */
  const Lorentz5Momentum& momentum() const { return theMomentum; }

  /**
   * Return the scale
   */
  Energy2 scale() const { return theScale; }

public:

  /**
   * Return a matrix element of the diagonal part, by contracting with
   * the given polarization vectors.
   */
  Complex matrixElement(const LorentzVector<Complex>& left,
			const LorentzVector<Complex>& right) const {
    return
      left.dot(momentum())*right.dot(momentum())/scale();
  }

};

}

#endif // HERWIG_SpinCorrelationTensor_H
