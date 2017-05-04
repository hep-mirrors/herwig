// -*- C++ -*-
//
// ShowerAlphaQED.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerAlphaQED_H
#define HERWIG_ShowerAlphaQED_H
//
// This is the declaration of the ShowerAlphaQED class.
//

#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value(scale) for \f$\alpha_{\rm QED}\f$.
 *  N.B. as we always use \f$\alpha(0)\f$ for the radiation of photons
 *  this class is very simple.
 *
 *  @see ShowerAlpha
 *  @see ShowerAlphaQCD
 *
 * @see \ref ShowerAlphaQEDInterfaces "The interfaces"
 * defined for ShowerAlphaQED.
 */
class ShowerAlphaQED: public ShowerAlpha {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerAlphaQED() : ShowerAlpha(), _alpha(1./137.) {}
  //@}

public:

  /**
   * Methods to return the coupling.
   * The methods are equivalent to the QCD ones 
   * and are necessary to make use of the virtuality of ShowerAlpha
   * at other places. 
   */
  //@{
  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   * @param scale The scale
   * @return The coupling
   */
  virtual double value(const Energy2 scale) const;

  /**
   * It returns the overestimiate of the coupling
   * multiplied by the scale factor scaleFactor().
   */
  virtual double overestimateValue() const;

  /**
   *  Return the ratio of the coupling at the scale to the overestimated value
   */
  virtual double ratio(const Energy2 scale,double factor=1.) const;
  //@}
  
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerAlphaQED> initShowerAlphaQED;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlphaQED & operator=(const ShowerAlphaQED &);

private:

  /**
   *  The value of the coupling, as we are producing real photons
   *  this is always \f$\alpha(q^2=0)\f$.
   */
  double _alpha; 
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerAlphaQED. */
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQED,1> {
  /** Typedef of the first base class of ShowerAlphaQED. */
  typedef Herwig::ShowerAlpha NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerAlphaQED class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerAlphaQED>
  : public ClassTraitsBase<Herwig::ShowerAlphaQED> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerAlphaQED"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerAlphaQED is implemented. It may also include several, space-separated,
   * libraries if the class ShowerAlphaQED depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ShowerAlphaQED_H */
