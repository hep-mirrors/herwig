// -*- C++ -*-
#ifndef HERWIG_ShowerAlpha_H
#define HERWIG_ShowerAlpha_H
//
// This is the declaration of the ShowerAlpha class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "ShowerVariables.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is the abstract class from which all types of running alphas 
 *  used in the Showering (Initial State alphaQCD, alphaQED, alphaEWK,...
 *  and Final State alphaQCD, alphaQED, alphaEWK,...) derived from.
 *  The main purpose of this class, and the ones that derive from, is
 *  to allow systematics evaluation of Initial State Radiation and,
 *  independently, of Final State Radiation effects, by allowing a
 *  multiplicative adimensional factor, which is 1.0 for the "central value",
 *  for the scale argument, mu^2, of the running alpha:
 *    <I> ShowerAlpha( mu^2 ) ---> ShowerAlpha( scaleFactor * mu^2 ) </I>
 *  This scale factor is a parameter which is settable by the user. 
 *  Although, of course, it is not clear my how much we should scale
 *  in order to get a "one-sigma" systematics error (but factors:
 *  1/2 and 2 are quite common), this method provides a double side error,
 *  and it appears more sensible than the rough and one-side evaluation
 *  via turning off the I.S.R. and/or F.S.R. (possibilities which are,
 *  anyway, provided by Herwig++). 
 *
 * @see ShowerAlphaQCD
 */ 
class SplittingGenerator; 

class ShowerAlpha: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  friend class SplittingGenerator;
  inline ShowerAlpha();
  inline ShowerAlpha(const ShowerAlpha &);
  virtual ~ShowerAlpha();

  /**
   * Pure virtual method that is supposed to return the 
   * running alpha value evaluated at the input scale.
   */
  virtual double value(const Energy2 scale) = 0;

  /**
   * Virtual method, which returns by default value( scale ), 
   * and that could be overrided in a derived class in the case an 
   * overestimate approximation of the alpha value is provided. 
   */
  virtual double overestimateValue() = 0;
  
  /**
   * It returns the adimensional factor that multiplies the 
   * scale argument, mu^2, of the running alpha:
   *    <I> ShowerAlpha( scaleFactor() * mu^2 ) </I>    
   * This is supposed to be <I>1.0</I> in "normal" conditions (central values) 
   * whereas different values can be useful for systematics evaluation 
   * for Initial State radiation or Final State radiation effects.
   */
  inline double scaleFactor() const;

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<ShowerAlpha> initShowerAlpha;

  /**
   * Private and non-existent assignment operator.
   */
  ShowerAlpha & operator=(const ShowerAlpha &);

  double _scaleFactor;

protected:

  ShowerVarsPtr _pointerShowerVariables;
  void setSV(ShowerVarsPtr scp); 

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerAlpha.
 */
template <>
struct BaseClassTrait<Herwig::ShowerAlpha,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ShowerAlpha>: public ClassTraitsBase<Herwig::ShowerAlpha> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ShowerAlpha"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "ShowerAlpha.icc"

#endif /* HERWIG_ShowerAlpha_H */
