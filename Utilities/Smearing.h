// -*- C++ -*-
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
     * Pure static class.
     */
    Smearing();
    Smearing(const Smearing & x);
    Smearing & operator=(const Smearing & x);

  };

}

#include "Smearing.icc"

#endif /* HERWIG_Smearing_H */

