// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQCD_H
#define HERWIG_ShowerAlphaQCD_H
//
// This is the declaration of the ShowerAlphaQCD class.
//

#include "ShowerAlpha.h"
#include "ShowerIndex.h"
#include "ShowerAlphaQCD.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value() and overestimateValue() for the 
 *  strong coupling.
 *
 *  @see ShowerAlpha
 *
 * @see \ref ShowerAlphaQCDInterfaces "The interfaces"
 * defined for ShowerAlphaQCD.
 */
class ShowerAlphaQCD: public ShowerAlpha {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerAlphaQCD();

  /**
   * The copy constructor.
   */
  inline ShowerAlphaQCD(const ShowerAlphaQCD &);

  /**
   * The destructor.
   */
  virtual ~ShowerAlphaQCD();
  //@}

public:

  /**
   *  Methods to return the coupling
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
  virtual double overestimateValue();
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
   *  Member function which calculate the coupling
   */
  //@{
  /**
   * The two-loop parametrization of \f$\alpha_S\f$.
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
  double alphaTwoLoop(Energy q, Energy lam, short nf); 

  /**
   * Return the two loop \f$\Lambda\f$ at the scale
   * Hacked in masses by hand for the moment before proper
   * interfacing...  obtained lambda solutions numerically in
   * Mathematica with my alphas.m using two-loop alphas from PDT2002
   * and as(M_Z=91.2GeV) = 0.118 *** ACHTUNG! *** this HAS to be done
   * automatically acc to the masses and as(M_Z) given by the PDT
   * class (which is supposed to be up-to-date).
   * @param q The scale
   * @return The number of flavours at the scale and \f$\Lambda\f$.
   */
  pair<short, Energy> getLamNfTwoLoop(Energy q); 

  /**
   * A toy parametrization of \f$\alpha_S\f$ with different parametrizations
   * of the IR behaviour, below q2min, set by type.  
   * Default is type = 1, i.e. alpha_s = 0 below q2min.
   * @param q2 The scale
   * @param q2min The minimum value below which the coupling is considered to 
   * be non-pertburative
   * @param type The type of behaviour
   */
  double alpha_s(double q2, double q2min, int type);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerAlphaQCD> initShowerAlphaQCD;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlphaQCD & operator=(const ShowerAlphaQCD &);

private:

  /**
   *  Minimum value of the scale
   */
  Energy _Qmin;

  /**
   *  Parameter controlling the behaviour of \f$\alpha_S\f$ in the non-perturbative
   *  region.
   */ 
  int _asType;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerAlphaQCD. */
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQCD,1> {
  /** Typedef of the first base class of ShowerAlphaQCD. */
  typedef Herwig::ShowerAlpha NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerAlphaQCD class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerAlphaQCD>
  : public ClassTraitsBase<Herwig::ShowerAlphaQCD> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ShowerAlphaQCD"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the ShowerAlphaQCD class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ShowerAlphaQCD.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQCD.tcc"
#endif

#endif /* HERWIG_ShowerAlphaQCD_H */
