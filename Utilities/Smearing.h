// -*- C++ -*-
//
// Smearing.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef HERWIG_Smearing_H
#define HERWIG_Smearing_H
//
// This is the declaration of the Smearing class.

#include "Herwig++/Config/Herwig.h"


namespace Herwig {

  using namespace ThePEG;

  /** \ingroup Utilities
   *
   * This is a pure static class which provides some useful methods for smearing.
   */
  class Smearing {

  public:

    /**
     * It returns true if succeed in drawing a value (x) drawn from a 
     * gaussian of specified mean and sigma, false otherwise. 
     * Indeed, it generates uncorrelated pairs and throws one of them away.
     */
    static bool gaussianSmearing( const double mean, const double sigma, double & x );

    /**
     * It returns true if it succeed in drawing a rotated 2-vector 
     * (vx,vy)  of given length  r , false otherwise.
     */
    static bool azimuthalSmearing( const double r, double & vx, double & vy );
    
  private:

    /**
     * Pure static class so no default constructor
     */
    Smearing();

    /**
     * Pure static class so no copy constructor
     */
    Smearing(const Smearing & x);

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    Smearing & operator=(const Smearing & x);

  };

}

#include "Smearing.icc"

#endif /* HERWIG_Smearing_H */

