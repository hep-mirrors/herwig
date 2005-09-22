// -*- C++ -*-
#ifndef HERWIG_SplitFun1to2_H
#define HERWIG_SplitFun1to2_H
//
// This is the declaration of the SplitFun1to2 class.

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerIndex.h"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This is an abstract class which defines the common interface <BR>
 *  for all <I>1-&GT;2</I> splitting functions, both for Initial State <BR>
 *  and Final State radiation. 
 *
 *  @see SplitFun
 *  @see SplitFun1to3
 *  @see QtoQGSplitFun
 *  @see QtoQGammaSplitFun
 *  @see GtoGGSplitFun
 *  @see GtoQQbarSplitFun
 */
class SplittingFunction: public Interfaced {

public:

  /**
   * Standard ctors and dtor.
   */
  inline SplittingFunction(ShowerIndex::InteractionType a) 
    : Interfaced(), _interactionType(a) {}
  inline SplittingFunction(const SplittingFunction &x) 
    : Interfaced(x), _interactionType(x._interactionType) {}
  virtual ~SplittingFunction() {}

  /**
   * It Specifies the interaction type (QCD,QED,EWK,...) of the vertex 
   * <I>A -&GT; B + C</I>, the PDG id for the emitter (<I>A</I>),
   * its mass, the ids of the two emission products (<I>B</I> and <I>C</I>).
   * and their masses.
   */

  /**
   * Pure virtual methods, that must be provided by the derived
   * classes, which return the exact values of the splitting function
   * <I>0 -&GT; 1 + 2</I> evaluated in terms of some combinations of:
   * z variable, phi azimuthal angle, helicities of the three particles.  
   * Notice that, to be more general, the returned value is complex, 
   * although for standard QCD it is indeed real.  
   */
  virtual double P(const double, const Energy2, const IdList &) = 0;

  /**
   * Virtual methods that could be overrided in a derived class in the
   * case it is known some approximate expressions for the exact
   * expressions of the splitting function (the ones provided by the
   * pure virtual methods above), and such that the former are always
   * above than the latter: overestimateXXX(...) | &GT;= |
   * XXX(...) |. The default definition of overestimateXXX(...) 
   * simply returns XXX(...).  Notice that we left the overestimated
   * functions and hence their integrals and inverses not depend on
   * qtilde2 since these are supposed to be as simple as possible.
   */
  virtual double overestimateP(const double z, const IdList &) = 0;  

  /**
   * The indefinite integral of the overestimated splitting function
   * overestimateIntegratedFun(z) and its inverse. 
   */
  virtual double integOverP(const double z) = 0; 
  virtual double invIntegOverP(const double r) = 0; 

  /**
   * The default definition of this method does nothing.
   * However, splitting functions which derive from this
   * class and that involve coloured particles should override 
   * this method, in order to provide the proper colour connection 
   * between the emitting parent and the branching products. 
   * More precisely, the first argument is a pair of pointers to 
   * ShowerColourLine objects, which are associated with, 
   * respectively, the colour (first element of the pair) and 
   * anticolour (second element of the pair) of the emitting particle.
   * The other two arguments of the method are also pair of pointers
   * to ShowerColourLine objects, for respectively the two 
   * branching products. Again, for each pair, the first element
   * is associated with the colour line and the second element
   * is associated with the anticolur line.
   * The ShowerColourLine objects pointed by any of the 
   * two elements of both pairs can be either one of object pointerd 
   * by the radiating parent, or can be a new one (and because of that, 
   * the pointers must be reference counted, not transient).
   * We prefer to use ShowerColourLine in this method
   * rather than ShowerParticles for two reason: first,
   * there is no coupling between SplittingFunction class
   * (and its derived classes) and the ShowerParticle class; 
   * and second, we avoid to repeat for each type of splitting
   * function the same operation of setting the color lines to
   * proper particles.
   */ 
  virtual void colourConnection(const ShoColinePair & parent,
				ShoColinePair & first,
				ShoColinePair & second) = 0;

  ShowerIndex::InteractionType interactionType() { return _interactionType; }

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  //static void Init(){}

private:

  static AbstractClassDescription<SplittingFunction> initSplittingFunction;
  SplittingFunction & operator=(const SplittingFunction &);
  
protected:

  ShowerIndex::InteractionType _interactionType;

};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.   
 */
template <>
struct BaseClassTrait<Herwig::SplittingFunction,1> {
  typedef ThePEG::Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SplittingFunction>
  : public ClassTraitsBase<Herwig::SplittingFunction>
{

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/SplittingFunction"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#endif /* HERWIG_SplittingFunction_H */
