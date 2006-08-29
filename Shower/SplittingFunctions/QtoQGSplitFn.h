// -*- C++ -*-
#ifndef HERWIG_QtoQGSplitFn_H
#define HERWIG_QtoQGSplitFn_H
//
// This is the declaration of the QtoQGSplitFn class.
//

#include "SplittingFunction.h"
#include "QtoQGSplitFn.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$q\to qg\f$. 
 *
 *  In this case the splitting function is given by
 * \f[P(z,t) =\frac{C_F}\left(\frac{1+z^2}{1-z}-2\frac{m^2_q}{t}\right),\f]
 * where \f$C_F=\frac43\f$.
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = \frac{2C_F}{1-z},\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z = -2C_F\ln(1-z),\f]
 * and its inverse is
 * \f[1-\exp\left(\frac{r}{2C_F}\right).\f]
 *
 * @see \ref QtoQGSplitFnInterfaces "The interfaces"
 * defined for QtoQGSplitFn.
 */
class QtoQGSplitFn: public SplittingFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline QtoQGSplitFn();
  //@}

public:

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
				ColinePair & second) const;

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
  static NoPIOClassDescription<QtoQGSplitFn> initQtoQGSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QtoQGSplitFn & operator=(const QtoQGSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QtoQGSplitFn. */
template <>
struct BaseClassTrait<Herwig::QtoQGSplitFn,1> {
  /** Typedef of the first base class of QtoQGSplitFn. */
  typedef Herwig::SplittingFunction NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QtoQGSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QtoQGSplitFn>
  : public ClassTraitsBase<Herwig::QtoQGSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::QtoQGSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QtoQGSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class QtoQGSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "QtoQGSplitFn.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoQGSplitFn.tcc"
#endif

#endif /* HERWIG_QtoQGSplitFn_H */
