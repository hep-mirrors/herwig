// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQED_H
#define HERWIG_ShowerAlphaQED_H
//
// This is the declaration of the <!id>ShowerAlphaQED<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This concrete class provides the definition of the 
// pure virtual function <!id>value(scale)<!!id>.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerAlpha.html">ShowerAlpha.h</a>,
// 

#include "ShowerAlpha.h"

namespace Herwig {

using namespace Pythia7;

class ShowerAlphaQED: public ShowerAlpha {

public:

  inline ShowerAlphaQED();
  inline ShowerAlphaQED(const ShowerAlphaQED &);
  virtual ~ShowerAlphaQED();
  // Standard ctors and dtor.

  // the following two methods are equivalent to the QCD ones 
  // and are necessary to make use of the virtuality of ShowerAlpha
  // at other places. 
  virtual double value(const Energy2 scale);
  // It returns the running alpha value evaluated at the input scale
  // multiplied by the scale factor <!id>scaleFactor()<!!id>.

  inline virtual double overestimateValue();
  // It returns the running alpha value evaluated at the input scale
  // multiplied by the scale factor <!id>scaleFactor()<!!id>.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

private:

  static ClassDescription<ShowerAlphaQED> initShowerAlphaQED;
  // Describe an concrete class with persistent data.

  ShowerAlphaQED & operator=(const ShowerAlphaQED &);
  //  Private and non-existent assignment operator.

private: 
  double _alpha; 
  
};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ShowerAlphaQED.
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQED,1> {
  typedef Herwig::ShowerAlpha NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ShowerAlphaQED>: public ClassTraitsBase<Herwig::ShowerAlphaQED> {
  static string className() { return "/Herwig++/ShowerAlphaQED"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ShowerAlphaQED.icc"

#endif /* HERWIG_ShowerAlphaQCD_H */
