// -*- C++ -*-
#ifndef HERWIG_SplittingFunction_H
#define HERWIG_SplittingFunction_H
//
// This is the declaration of the SplittingFunction class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower2/ShowerConfig.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "Herwig++/Shower2/Couplings/ShowerIndex.h"
#include "SplittingFunction.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This is an abstract class which defines the common interface
 *  for all \f$1\to2\f$ splitting functions, for both initial-state
 *  and final-state radiation. 
 *
 *  The SplittingFunction class contains a number of purely virtual members
 *  which must be implemented in the inheriting classes. The class also stores
 *  the interaction type of the spltting function.
 *
 *  The inheriting classes need to specific the splitting function 
 *  \f$P(z,\tilde{q}^2)\f$, in terms of the energy fraction \f$z\f$ and evolution scale
 *  \f$\tilde{q}^2\f$ for a given splitting. In addition an overestimate of the 
 *  splitting function, \f$P_{\rm over}(z)\f$ which only depends upon \f$z\f$, 
 *  the integral, inverse of the integral for this overestimate and
 *  ratio of the true splitting function to the overestimate must be provided
 *  as they are necessary for the veto alogrithm used to implement the evolution.
 *
 * @see \ref SplittingFunctionInterfaces "The interfaces"
 * defined for SplittingFunction.
 */
class SplittingFunction: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.   
   * @param a All splitting functions must have an interaction type
   * @param b All splitting functions must have an interaction order
   */
  inline SplittingFunction(ShowerIndex::InteractionType a, unsigned int b);

  /**
   * The copy constructor.
   */
  inline SplittingFunction(const SplittingFunction &);

  /**
   * The destructor.
   */
  virtual ~SplittingFunction();
  //@}

public:

  /**
   *  Methods to return the interaction type and order for the splitting function
   */
  //@{
  /**
   *  Return the type of the interaction
   */
  inline ShowerIndex::InteractionType interactionType();

  /**
   *  Return the order of the splitting function in the interaction
   */
  inline unsigned int interactionOrder();
  //@}

  /**
   *  Purely virtual method which should determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) = 0;

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * Purely virtual method which should return the exact value of the splitting function,
   * \f$P\f$ evaluated in terms of the energy fraction, \f$z\f$, and the evolution scale 
   \f$\tilde{q}^2\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids) = 0;

  /**
   * Purely virtual method which should return
   * an overestimate of the splitting function,
   * \f$P_{\rm over}\f$ such that the result \f$P_{\rm over}\geq P\f$. This function
   * should be simple enough that it does not depend on the evolution scale.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList & ids) = 0; 

  /**
   * Purely virtual method which should return
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids) = 0;

  /**
   * Purely virtual method which should return the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   */
  virtual double integOverP(const double z) = 0; 

  /**
   * Purely virtual method which should return the inverse of the 
   * indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$ which is used to
   * generate the value of \f$z\f$.
   * @param r Value of the splitting function to be inverted
   */
  virtual double invIntegOverP(const double r) = 0; 
  //@}

  /**
   * Purely virtual method which should make the proper colour connection 
   * between the emitting parent and the branching products.
   *
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
   *
   * The ColourLine objects pointed by any of the 
   * two elements of both pairs can be either one of object pointerd 
   * by the radiating parent, or can be a new one (and because of that, 
   * the pointers must be reference counted, not transient).
   * We prefer to use ColourLine in this method
   * rather than ShowerParticles for two reasons: first,
   * there is no coupling between SplittingFunction class
   * (and its derived classes) and the ShowerParticle class; 
   * and second, we avoid repeating for each type of splitting
   * function the same operation of setting the color lines to
   * proper particles.
   */ 
  virtual void colourConnection(const ColinePair & parent,
				ColinePair & first,
				ColinePair & second) = 0;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<SplittingFunction> initSplittingFunction;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SplittingFunction & operator=(const SplittingFunction &);

private:

  /**
   *  The interaction type for the splitting function.
   */
  ShowerIndex::InteractionType _interactionType;

  /**
   *  The order of the splitting function in the coupling
   */
  unsigned int _interactionorder;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SplittingFunction. */
template <>
struct BaseClassTrait<Herwig::SplittingFunction,1> {
  /** Typedef of the first base class of SplittingFunction. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SplittingFunction class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SplittingFunction>
  : public ClassTraitsBase<Herwig::SplittingFunction> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SplittingFunction"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SplittingFunction is implemented. It may also include several, space-separated,
   * libraries if the class SplittingFunction depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "SplittingFunction.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SplittingFunction.tcc"
#endif

#endif /* HERWIG_SplittingFunction_H */
