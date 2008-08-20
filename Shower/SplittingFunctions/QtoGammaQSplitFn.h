// -*- C++ -*-
#ifndef HERWIG_QtoGammaQSplitFn_H
#define HERWIG_QtoGammaQSplitFn_H
//
// This is the declaration of the QtoGammaQSplitFn class.
//

#include "SplittingFunction.h"
#include "QtoGammaQSplitFn.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The QtoGammaQSplitFn class
 * provides the concrete implementation of the exact leading-order 
 * splitting function for \f$q\to \gamma q\f$.
 *
 *  In this case the splitting function is given by
 * \f[P(z,t) = \left(\frac{2(1-z)+z^2}{z}-2\frac{m^2_q}t\right).\f]
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = 2\frac1z,\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z = 2\ln z,\f]
 * and its inverse is
 * \f[\exp\left(\frac{r}{2}\right).\f]
 *
 *  @see SplittingFunction
 *
 */
class QtoGammaQSplitFn: public SplittingFunction {

public:

  /**
   * The default constructor.
   */
  inline QtoGammaQSplitFn();


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
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double integOverP(const double z, const IdList & ids, 
			       unsigned int PDFfactor=0) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r, const IdList & ids, 
			       unsigned int PDFfactor=0) const;
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
  static NoPIOClassDescription<QtoGammaQSplitFn> initQtoGammaQSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QtoGammaQSplitFn & operator=(const QtoGammaQSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QtoGammaQSplitFn. */
template <>
struct BaseClassTrait<Herwig::QtoGammaQSplitFn,1> {
  /** Typedef of the first base class of QtoGammaQSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QtoGammaQSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QtoGammaQSplitFn>
  : public ClassTraitsBase<Herwig::QtoGammaQSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QtoGammaQSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QtoGammaQSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class QtoGammaQSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "QtoGammaQSplitFn.icc"

#endif /* HERWIG_QtoGammaQSplitFn_H */
