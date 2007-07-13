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

/** \ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$g\to gg\f$. 
 *
 * In this case the splitting function is given by
 * \f[P(z) =2C_A*\left(\frac{z}{1-z}+\frac{1-z}{z}+z(1-z)\right),\f]
 * where \f$C_A=3\f$
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = 2C_A\left(\frac1z+\frac1{1-z}\right),\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z =2C_a\ln\left(\frac{z}{1-z}\right),\f]
 * and its inverse is
 * \f[\frac{\exp\left(\frac{r}{2C_A}\right)}{(1+\exp\left(\frac{r}{2C_A}\right)}\f]
 *
 *
 * @see \ref GtoGGSplitFnInterfaces "The interfaces"
 * defined for GtoGGSplitFn.
 */
class GtoGGSplitFn: public SplittingFunction {

public:

  /**
   * The default constructor.
   */
  inline GtoGGSplitFn();

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const;

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * The concrete implementation of the splitting function, \f$P(z,t)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids,
		   const bool mass) const;

  /**
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList & ids) const; 

  /**
   * The concrete implementation of the
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,t)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			const bool mass) const;

  /**
   * The concrete implementation of the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   */
  virtual double integOverP(const double z) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   */ 
  virtual double invIntegOverP(const double r) const;
  //@}

  /**
   * Purely virtual method which should make the proper colour connection 
   * between the emitting parent and the branching products.
   * @param parent The parent for the branching
   * @param first  The first  branching product
   * @param second The second branching product
   * @param back Whether this is foward or backward evolution.
   */
  virtual void colourConnection(tShowerParticlePtr parent,
				tShowerParticlePtr first,
				tShowerParticlePtr second,
				const bool back) const;

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
  static string className() { return "Herwig::GtoGGSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GtoGGSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class GtoGGSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "GtoGGSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GtoGGSplitFn.tcc"
#endif

#endif /* HERWIG_GtoGGSplitFn_H */
