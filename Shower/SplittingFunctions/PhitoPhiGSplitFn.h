// -*- C++ -*-
//
// PhitoPhiGSplitFn.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PhitoPhiGSplitFn_H
#define HERWIG_PhitoPhiGSplitFn_H
//
// This is the declaration of the PhitoPhiGSplitFn class.
//

#include "SplittingFunction.h"
#include "PhitoPhiGSplitFn.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$\phi\to \phi g\f$.
 *
 * In this case the splitting function is given by
 * \f[P(z,t) = 2C_F\left(\frac{z}{1-z}-\frac{m^2_\phi}{t}\right),\f]
 * where \f$C_F=\frac43\f$.
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = \frac{2C_F}{1-z},\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z = -2C_F\ln(1-z),\f]
 * and its inverse is
 * \f[1-\exp\left(\frac{r}{2C_F}\right).\f]
 *
 * @see \ref PhitoPhiGSplitFnInterfaces "The interfaces"
 * defined for PhitoPhiGSplitFn.
 */
class PhitoPhiGSplitFn: public SplittingFunction {

public:

  /**
   * The default constructor.
   */
  inline PhitoPhiGSplitFn();

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
   * The concrete implementation of the splitting function, \f$P\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids,
		   bool mass) const;

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
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			bool mass) const;

  /**
   * The concrete implementation of the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double integOverP(const double z, unsigned int PDFfactor=0) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r, unsigned int PDFfactor=0) const;
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
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<PhitoPhiGSplitFn> initPhitoPhiGSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PhitoPhiGSplitFn & operator=(const PhitoPhiGSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PhitoPhiGSplitFn. */
template <>
struct BaseClassTrait<Herwig::PhitoPhiGSplitFn,1> {
  /** Typedef of the first base class of PhitoPhiGSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PhitoPhiGSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PhitoPhiGSplitFn>
  : public ClassTraitsBase<Herwig::PhitoPhiGSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PhitoPhiGSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PhitoPhiGSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class PhitoPhiGSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "PhitoPhiGSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PhitoPhiGSplitFn.tcc"
#endif

#endif /* HERWIG_PhitoPhiGSplitFn_H */
