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
 * \f[P(z,\tilde{q}^2) = \frac{C_F}{z}\left( 2(1-z)+z^2-2\frac{m^2_q}{\tilde{q}^2(1-z)}
 * \right),\f]
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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline QtoGQSplitFn();

  /**
   * The copy constructor.
   */
  inline QtoGQSplitFn(const QtoGQSplitFn &);

  /**
   * The destructor.
   */
  virtual ~QtoGQSplitFn();
  //@}

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids);

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
   * The concrete implementation of the
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids);

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
   * @param parent Pair of pointers to ColourLine objects, 
   * which are associated with, 
   * respectively, the colour (first element of the pair) and 
   * anticolour (second element of the pair) of the emitting particle.
   * @param first Pair of pointers
   * to ColourLine objects, for respectively the first 
   * branching product. Again the first element
   * is associated with the colour line and the second element
   * is associated with the anticolur line.
   * @param second As first but for the second particle.
   */
  virtual void colourConnection(const ColinePair & parent,
				ColinePair & first,
				ColinePair & second);

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
  static string className() { return "Herwig++::QtoGQSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QtoGQSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class QtoGQSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "QtoGQSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoGQSplitFn.tcc"
#endif

#endif /* HERWIG_QtoGQSplitFn_H */
