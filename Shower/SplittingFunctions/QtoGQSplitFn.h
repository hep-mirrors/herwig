// -*- C++ -*-
#ifndef HERWIG_QtoGQSplitFn_H
#define HERWIG_QtoGQSplitFn_H
//
// This is the declaration of the QtoGQSplitFn class.
//

#include "SplittingFunction.h"
#include "QtoGQSplitFn.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This classs provides the concrete implementation of the exact leading-order 
 *  splitting function for \f$q\to gq\f$.
 *
 *  In this case the splitting function is given by
 * \f[P(z,t) = C_F\left(\frac{2(1-z)+z^2}{z}-2\frac{m^2_q}t\right),\f]
 * where \f$C_F=\frac43\f$.
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = 2C_F\frac1z,\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z = 2C_F\ln z,\f]
 * and its inverse is
 * \f[\exp\left(\frac{r}{2C_F}\right).\f]
 *
 *  @see SplittingFunction
 */
class QtoGQSplitFn: public SplittingFunction {

public:

  /**
   * The default constructor.
   */
  inline QtoGQSplitFn();

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

  /**
   * Purely virtual method which should return the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}/(1-z)\f$.
   * @param z   The energy fraction.
   */
  virtual double integOverPPDFFactor(const double z) const; 

  /**
   * Purely virtual method which should return the inverse of the 
   * indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}/(1-z)\f$ which is used to
   * generate the value of \f$z\f$.
   * @param r Value of the splitting function to be inverted
   */
  virtual double invIntegOverPPDFFactor(const double r) const; 
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
  static NoPIOClassDescription<QtoGQSplitFn> initQtoGQSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QtoGQSplitFn & operator=(const QtoGQSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QtoGQSplitFn. */
template <>
struct BaseClassTrait<Herwig::QtoGQSplitFn,1> {
  /** Typedef of the first base class of QtoGQSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QtoGQSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QtoGQSplitFn>
  : public ClassTraitsBase<Herwig::QtoGQSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QtoGQSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QtoGQSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class QtoGQSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "QtoGQSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoGQSplitFn.tcc"
#endif

#endif /* HERWIG_QtoGQSplitFn_H */
