// -*- C++ -*-
//
#ifndef HERWIG_Smearing_H
#define HERWIG_Smearing_H
//
// This is the declaration of the <!id>Smearing<!!id> class.
//
// CLASSDOC SUBSECTION Description:
// 
// This is a pure static class which provides some useful methods for smearing.
//

#include "Herwig++/Config/Herwig.h"


namespace Herwig {


  using namespace Pythia7;


  class Smearing {

  public:

    static bool gaussianSmearing( const double mean, const double sigma, double & x );
    // It returns true if succeed in drawing a value (x) drawn from a 
    // gaussian of specified <!id>mean<!!id> and <!id>sigma<!!id>,
    // false otherwise. 
    // Indeed, it generates uncorrelated pairs and throws one of them away.

    static bool azimuthalSmearing( const double r, double & vx, double & vy );
    // It returns true if it succeed in drawing a rotated 2-vector 
    // <!id>(vx,vy)<!!id> of given length <!id>r<!!id>, false otherwise.
    
  private:

    Smearing();
    Smearing(const Smearing & x);
    Smearing & operator=(const Smearing & x);
    // pure static class

  };


}

#include "Smearing.icc"

#endif /* HERWIG_Smearing_H */

