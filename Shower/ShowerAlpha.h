// -*- C++ -*-
#ifndef HERWIG_ShowerAlpha_H
#define HERWIG_ShowerAlpha_H
//
// This is the declaration of the <!id>ShowerAlpha<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the abstract class from which all types of running alphas <BR>
// used in the Showering (Initial State alphaQCD, alphaQED, alphaEWK,... <BR>
// and Final State alphaQCD, alphaQED, alphaEWK,...) derived from. <BR>
// The main purpose of this class, and the ones that derive from, is <BR>
// to allow systematics evaluation of Initial State Radiation and, <BR>
// independently, of Final State Radiation effects, by allowing a <BR>
// multiplicative adimensional factor, which is 1.0 for the "central value", <BR>
// for the scale argument, mu^2, of the running alpha: <BR>
//   <I> ShowerAlpha( mu^2 ) ---> ShowerAlpha( scaleFactor * mu^2 ) </I><BR>
// This scale factor is a parameter which is settable by the user. <BR>
// Although, of course, it is not clear my how much we should scale <BR>
// in order to get a "one-sigma" systematics error (but factors: <BR>
// 1/2 and 2 are quite common), this method provides a double side error, <BR>
// and it appears more sensible than the rough and one-side evaluation <BR>
// via turning off the I.S.R. and/or F.S.R. (possibilities which are, <BR>
// anyway, provided by Herwig++). 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerAlphaQCD.html">ShowerAlphaQCD.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "ShowerConfig.h"

namespace Herwig {

using namespace Pythia7;

class ShowerAlpha: public Pythia7::HandlerBase {

public:

  inline ShowerAlpha();
  inline ShowerAlpha(const ShowerAlpha &);
  virtual ~ShowerAlpha();
  // Standard ctors and dtor.

  virtual double value(const Energy2 scale) = 0;
  // Pure virtual method that is supposed to return the 
  // running alpha value evaluated at the input scale.

  virtual double overestimateValue() = 0;
  // Virtual method, which returns by default value( scale ), 
  // and that could be overrided in a derived class in the case an 
  // overestimate approximation of the alpha value is provided. 
  
  inline double scaleFactor() const;
  // It returns the adimensional factor that multiplies the 
  // scale argument, mu^2, of the running alpha:
  //   <I> ShowerAlpha( scaleFactor() * mu^2 ) </I>    
  // This is supposed to be <I>1.0</I> in "normal" conditions (central values) 
  // whereas different values can be useful for systematics evaluation 
  // for Initial State radiation or Final State radiation effects.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static AbstractClassDescription<ShowerAlpha> initShowerAlpha;
  // Describe an abstract base class with persistent data.

  ShowerAlpha & operator=(const ShowerAlpha &);
  //  Private and non-existent assignment operator.

  double _scaleFactor;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ShowerAlpha.
template <>
struct BaseClassTrait<Herwig::ShowerAlpha,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ShowerAlpha>: public ClassTraitsBase<Herwig::ShowerAlpha> {
  static string className() { return "/Herwig++/ShowerAlpha"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ShowerAlpha.icc"

#endif /* HERWIG_ShowerAlpha_H */
