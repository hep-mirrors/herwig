// -*- C++ -*-
//
// ShowerAlpha.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerAlpha_H
#define HERWIG_ShowerAlpha_H
//
// This is the declaration of the ShowerAlpha class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerAlpha.fh"

namespace Herwig {

using namespace ThePEG;

/**  \ingroup Shower
 * 
 *  This class is the abstract class from which all types of running couplings 
 *  used in the Showering derive from.
 *  The main purpose of this class, and the ones that derive from it, is
 *  to allow systematic uncertainties for the initial-state radiation and,
 *  independently, the final-state radiation effects, to be evaluated.
 *
 *  This is achieved by allowing a multiplicative factor,
 *  which is 1.0 for the "central value",
 *  for the scale argument, \f$\mu^2\f$, of the running coupling. This
 *  factor, \f$f\f$ is given by the scaleFactor() member and the coupling
 *  returned by the value() member is such that
 * \f[\alpha(\mu^2)\to \alpha(f\times\mu^2).\f]
 *  This scale factor is a parameter which is settable by the user, via the
 * interface. 
 *  Although, of course, it is not clear my how much we should scale
 *  in order to get a \f$1\sigma\f$ systematic error (but factors:
 *  1/2 and 2 are quite common), this method provides a double side error,
 *  and it appears more sensible than the rough and one-sided evaluation
 *  obtained
 *  via turning off the I.S.R. and/or F.S.R. (possibilities which are,
 *  anyway, provided by Herwig++). 
 *
 * @see \ref ShowerAlphaInterfaces "The interfaces"
 * defined for ShowerAlpha.
 */
class ShowerAlpha: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerAlpha();
  //@}

public:

  /**
   *  Methods to return the coupling and the scaleFactor
   */
  //@{
  /**
   * Pure virtual method that is supposed to return the 
   * running alpha value evaluated at the input scale.
   * @param scale The scale
   * @return The coupling
   */
  virtual double value(const Energy2 scale) const = 0;

  /**
   * Virtual method, which 
   * should be overridden in a derived class to provide an 
   * overestimate approximation of the alpha value. 
   */
  virtual double overestimateValue() const = 0;

  /**
   *  Virtual method which returns the ratio of the running alpha
   * value at the input scale to the overestimated value.
   * @param scale The scale
   * @return The ratio
   */
  virtual double ratio(const Energy2 scale) const = 0;

  /**
   * It returns the factor that multiplies the 
   * scale argument, \f$\mu^2\f$, of the running coupling.
   * This is supposed to be <I>1.0</I> in normal conditions (central values) 
   * whereas different values can be useful for systematics evaluation 
   * for Initial State radiation or Final State radiation effects.
   */
  inline double scaleFactor() const;

  /**
   * Initialize this coupling.
   */
  virtual inline void initialize ();

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<ShowerAlpha> initShowerAlpha;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlpha & operator=(const ShowerAlpha &);

private:

  /**
   *  The scale factor
   */
  double _scaleFactor;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerAlpha. */
template <>
struct BaseClassTrait<Herwig::ShowerAlpha,1> {
  /** Typedef of the first base class of ShowerAlpha. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerAlpha class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerAlpha>
  : public ClassTraitsBase<Herwig::ShowerAlpha> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerAlpha"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerAlpha is implemented. It may also include several, space-separated,
   * libraries if the class ShowerAlpha depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "ShowerAlpha.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlpha.tcc"
#endif

#endif /* HERWIG_ShowerAlpha_H */
