// -*- C++ -*-
//
// OneHalfHalfSplitFn.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_OneHalfHalfSplitFn_H
#define HERWIG_OneHalfHalfSplitFn_H
//
// This is the declaration of the OneHalfHalfSplitFn class.
//

#include "SplittingFunction.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$1\to \frac12\frac12\f$. 
 *
 * In this case the splitting function is given by
 * \f[P(z,t) =C\left(1-2z(1-z)+2\frac{m_q^2}{t}\right),\f]
 * where \f$C\f$ is the corresponding colour factor
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = C,\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z =Cz,\f]
 * and its inverse is
 * \f[\frac{r}{C}\f]
 *
 * @see \ref OneHalfHalfSplitFnInterfaces "The interfaces"
 * defined for OneHalfHalfSplitFn.
 */
class OneHalfHalfSplitFn: public SplittingFunction {

public:

  /**
   * The default constructor.
   */
  inline OneHalfHalfSplitFn() : SplittingFunction(1) {}

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
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
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
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double integOverP(const double z,  const IdList & ids, 
			    unsigned int PDFfactor=0) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r,  const IdList & ids, 
			       unsigned int PDFfactor=0) const;
  //@}

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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<OneHalfHalfSplitFn> initOneHalfHalfSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OneHalfHalfSplitFn & operator=(const OneHalfHalfSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of OneHalfHalfSplitFn. */
template <>
struct BaseClassTrait<Herwig::OneHalfHalfSplitFn,1> {
  /** Typedef of the first base class of OneHalfHalfSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the OneHalfHalfSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::OneHalfHalfSplitFn>
  : public ClassTraitsBase<Herwig::OneHalfHalfSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::OneHalfHalfSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * OneHalfHalfSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class OneHalfHalfSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_OneHalfHalfSplitFn_H */
