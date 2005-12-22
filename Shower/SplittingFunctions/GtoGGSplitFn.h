// -*- C++ -*-
#ifndef HERWIG_GtoGGSplitFn_H
#define HERWIG_GtoGGSplitFn_H
//
// This is the declaration of the GtoGGSplitFn class.
//

#include "SplittingFunction.h"
#include "GtoGGSplitFn.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$g\to gg\f$. 
 *
 *  In this case the splitting function is given by
 * \f[P(z,\tilde{q}^2) =2C_A*\left(\frac{z}{1-z}+\frac{1-z}{z}+z(1-z)\right),\f]
 * where \f$C_A=3\f$
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = 2C_A\left(\frac1z+\frac1{1-z}\right),\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z =2C_a\ln\left(\frac{z}{1-z}\right),\f]
 * and its inverse is
 * \f[\frac{\exp\left(\frac{r}{2C_A}\right)}{(1+\exp\left(\frac{r}{2C_A}\right)}\f]
 *
 *  @see SplittingFunction
 */
class GtoGGSplitFn: public SplittingFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline GtoGGSplitFn();

  /**
   * The copy constructor.
   */
  inline GtoGGSplitFn(const GtoGGSplitFn &);

  /**
   * The destructor.
   */
  virtual ~GtoGGSplitFn();
  //@}

public:

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * The concrete implementation of the splitting function, \f$P\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids);


  /**
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList & ids); 

  /**
   * The concrete implementation of the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   */
  virtual double integOverP(const double z);

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   */ 
  virtual double invIntegOverP(const double r);
  //@}

  /**
   *  Concrete implementation of the method to make the colour connections.
   * @param parent Pair of pointers to ShowerColourLine objects, 
   * which are associated with, 
   * respectively, the colour (first element of the pair) and 
   * anticolour (second element of the pair) of the emitting particle.
   * @param first Pair of pointers
   * to ShowerColourLine objects, for respectively the first 
   * branching product. Again the first element
   * is associated with the colour line and the second element
   * is associated with the anticolur line.
   * @param second As first but for the second particle.
   */
  virtual void colourConnection(const ShoColinePair & parent,
				ShoColinePair & first,
				ShoColinePair & second);

public:

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<GtoGGSplitFn> initGtoGGSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GtoGGSplitFn & operator=(const GtoGGSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GtoGGSplitFn. */
template <>
struct BaseClassTrait<Herwig::GtoGGSplitFn,1> {
  /** Typedef of the first base class of GtoGGSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GtoGGSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GtoGGSplitFn>
  : public ClassTraitsBase<Herwig::GtoGGSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GtoGGSplitFn"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GtoGGSplitFn class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "GtoGGSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GtoGGSplitFn.tcc"
#endif

#endif /* HERWIG_GtoGGSplitFn_H */
