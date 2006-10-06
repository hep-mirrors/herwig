// -*- C++ -*-
#ifndef HERWIG_Interpolator_H
#define HERWIG_Interpolator_H
//
// This is the declaration of the Interpolator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Interpolator.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Utilities
 *  \author Peter Richardson
 *
 *  This class implments a polynominal interpolation of a table of values, it is
 *  based on the interpolation code in FORTRAN HERWIG. 
 *
 * @see \ref InterpolatorInterfaces "The interfaces"
 * defined for Interpolator.
 */
class Interpolator: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Interpolator();

  /**
   * Constructor with data as vectors.
   */
  Interpolator(vector<double> f, vector<double> x, int order);
  //@}

  /**
   *  Return the interpolated value
   */
  double operator () (double) const;

public:

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Interpolator> initInterpolator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Interpolator & operator=(const Interpolator &);
  
private:
  
  /**
   * The x values.
   */
  vector<double> _xval;

  /**
   * the function values.
   */
  vector<double> _fun;

  /**
   * the order of interpolation.
   */
  unsigned int _order;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Interpolator. */
template <>
struct BaseClassTrait<Herwig::Interpolator,1> {
  /** Typedef of the first base class of Interpolator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Interpolator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Interpolator>
  : public ClassTraitsBase<Herwig::Interpolator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::Interpolator"; }
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

#include "Interpolator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Interpolator.tcc"
#endif

#endif /* HERWIG_Interpolator_H */
