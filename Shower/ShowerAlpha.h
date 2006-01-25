// -*- C++ -*-
#ifndef HERWIG_ShowerAlpha_H
#define HERWIG_ShowerAlpha_H
//
// This is the declaration of the ShowerAlpha class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerConfig.h"
#include "ShowerVariables.h"
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
 * interfaced. 
 *  Although, of course, it is not clear my how much we should scale
 *  in order to get a \f$1\sigma\f$ systematic error (but factors:
 *  1/2 and 2 are quite common), this method provides a double side error,
 *  and it appears more sensible than the rough and one-sided evaluation
 *  obtained
 *  via turning off the I.S.R. and/or F.S.R. (possibilities which are,
 *  anyway, provided by Herwig++). 
 *
 * @see ShowerAlphaQCD
 * @see ShowerAlphaQED
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

  /**
   * The copy constructor.
   */
  inline ShowerAlpha(const ShowerAlpha &);

  /**
   * The destructor.
   */
  virtual ~ShowerAlpha();
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
  virtual double value(const Energy2 scale) = 0;

  /**
   * Virtual method, which returns by default value( scale ), 
   * and that could be overrided in a derived class in the case an 
   * overestimate approximation of the alpha value is provided. 
   */
  virtual double overestimateValue() = 0;

  /**
   * It returns the factor that multiplies the 
   * scale argument, \f$\mu^2\f$, of the running coupling.
   * This is supposed to be <I>1.0</I> in normal conditions (central values) 
   * whereas different values can be useful for systematics evaluation 
   * for Initial State radiation or Final State radiation effects.
   */
  inline double scaleFactor() const;
  //@}

  /**
   *  Methods to set/get to the ShowerVariables
   */
  //@{
  /**
   *  Set the ShowerVariables pointer
   */
  inline void setSV(ShowerVarsPtr);

  /**
   *  Get the ShowerVariables pointer
   */
  inline tcShowerVarsPtr getSV();
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

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

  /**
   *  Pointer to the ShowerVariables
   */
  ShowerVarsPtr _pointerShowerVariables;
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
  static string className() { return "Herwig++::ShowerAlpha"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the ShowerAlpha class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ShowerAlpha.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlpha.tcc"
#endif

#endif /* HERWIG_ShowerAlpha_H */
