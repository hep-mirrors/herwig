// -*- C++ -*-
#ifndef HERWIG_SplittingFunction_H
#define HERWIG_SplittingFunction_H
//
// This is the declaration of the SplittingFunction class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "Herwig++/Shower/ShowerIndex.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "SplittingFunction.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This is an abstract class which defines the common interface
 *  for all \f$1\to2\f$ splitting functions, both for initial-state
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
 *  the integral and inverse of the integral for this overestimate must be provided
 *  as they are necessary for the veto alogrithm used to implement the evolution.
 *
 *  @see QtoGQSplitFn
 *  @see QtoQGSplitFn
 *  @see QtoQGammaSplitFn
 *  @see GtoGGSplitFn
 *  @see GtoQQbarSplitFn
 */
class SplittingFunction: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   * @param a All splitting functions must have an interaction type
   */
  inline SplittingFunction(ShowerIndex::InteractionType a);

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
   *  Return the type of the interaction
   */
  inline ShowerIndex::InteractionType interactionType();

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
   *
   * The ShowerColourLine objects pointed by any of the 
   * two elements of both pairs can be either one of object pointerd 
   * by the radiating parent, or can be a new one (and because of that, 
   * the pointers must be reference counted, not transient).
   * We prefer to use ShowerColourLine in this method
   * rather than ShowerParticles for two reasons: first,
   * there is no coupling between SplittingFunction class
   * (and its derived classes) and the ShowerParticle class; 
   * and second, we avoid repeating for each type of splitting
   * function the same operation of setting the color lines to
   * proper particles.
   */ 
  virtual void colourConnection(const ShoColinePair & parent,
				ShoColinePair & first,
				ShoColinePair & second) = 0;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

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
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SplittingFunction class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "SplittingFunction.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SplittingFunction.tcc"
#endif

#endif /* HERWIG_SplittingFunction_H */
