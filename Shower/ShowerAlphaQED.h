// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQED_H
#define HERWIG_ShowerAlphaQED_H
//
// This is the declaration of the ShowerAlphaQED class.
//

#include "ShowerAlpha.h"
#include "ShowerAlphaQED.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value(scale) for \f$\alpha_{\rm QED}\f$.
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
  inline ShowerAlphaQED();

  /**
   * The copy constructor.
   */
  inline ShowerAlphaQED(const ShowerAlphaQED &);

  /**
   * The destructor.
   */
  virtual ~ShowerAlphaQED();
  //@}

public:

  /**
   *  Methods to return the coupling.
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
  virtual double value(const Energy2 scale);

  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  inline virtual double overestimateValue();
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

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
  static string className() { return "Herwig++::ShowerAlphaQED"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the ShowerAlphaQED class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ShowerAlphaQED.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQED.tcc"
#endif

#endif /* HERWIG_ShowerAlphaQED_H */
