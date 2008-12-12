// -*- C++ -*-
//
// Interpolator2d.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Interpolator2d_H
#define HERWIG_Interpolator2d_H
//
// This is the declaration of the Interpolator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Config/Pointers.h"
#include <cassert>
#include <vector>
#include <math.h>
#include <iostream>

namespace Herwig {

  using namespace std;
  using namespace ThePEG;

  class Interpolator2d: public Interfaced {

  public:

    //Constructor of 2d interpolator object
    //takes a vector of the function weights evaluated on a m*n grid at grid points (i, j)
    //also takes vectors of the x1_i and x2_j values on that grid
    Interpolator2d( const vector< vector< double > > &, 
		    const vector< double > &,
		    const vector< double > & );

    /**
     *  Return the interpolated value
     */
    double operator () (double, double) const;
    

    /** @name Functions used by the persistent I/O system. */
    //@{
    /**
     * Function used to write out object persistently.
     * @param os the persistent output stream written to.
     */
    void persistentOutput(PersistentOStream & os) const;

    /**
     * Function used to read in object persistently.
     * @param is the persistent input stream read from.
     * @param version the version number of the object when written.
     */
    void persistentInput(PersistentIStream & is, int version);
    //@}

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

  protected:

    /** @name Clone Methods. */
    //@{
    /**
     * Make a simple clone of this object.
     * @return a pointer to the new object.
     */
    virtual IBPtr clone() const { return new_ptr(*this); }

    /** Make a clone of this object, possibly modifying the cloned object
     * to make it sane.
     * @return a pointer to the new object.
     */
    virtual IBPtr fullclone() const { return new_ptr(*this); }
    //@}

  private:

    //returns the 4*4 coefficient matrix required for grid squar i,j
    vector< vector< double > > getCoefficient( unsigned int, unsigned int );
    
    //functions to approximate the various derivatives required using the
    //centered difference approximation
    vector< vector< double > > findDerivX1();
    vector< vector< double > > findDerivX2();
    vector< vector< double > > findDerivX1X2();

    void gridHists();
    
    //the coefficient matrices - a 4x4 matrix for each grid square
    vector< vector< vector< vector< double > > > > _coefficients;
    
    //the (approximate) partitial derivatives at each grid point
    vector< vector< double > > _theDerivX1;
    vector< vector< double > > _theDerivX2;
    vector< vector< double > > _theDerivX1X2;
    
    //function weights on grid x1_i x2_j
    vector< vector< double > > _the_weights;
    //grid x1_i points
    vector< double > _x1_points;
    //grid x2_j points
    vector< double > _x2_points;
    //grid separation (must be constant) in x1 and x2 directions
    double  _dx1, _dx2;
    //magic matrix from numerical recipes
    vector< vector< int > > _wt;

    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is a concrete class with persistent data.
     */
    static ClassDescription<Interpolator2d > initInterpolator;  
  };
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  ThePEG_DECLARE_POINTERS(Herwig::Interpolator2d, Interpolator2dPtr);


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Interpolator. */
template <>
struct BaseClassTrait< Herwig::Interpolator2d, 1 > {
  /** Typedef of the first base class of Interpolator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Interpolator class and the shared object where it is defined. */
template <>
struct ClassTraits< Herwig::Interpolator2d >
  : public ClassTraitsBase< Herwig::Interpolator2d > {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Interpolator2d"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Interpolator is implemented. It may also include several, space-separated,
   * libraries if the class Interpolator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwUtils.so"; }
};

/** @endcond */

}

#endif /* HERWIG_Interpolator2d_H */
