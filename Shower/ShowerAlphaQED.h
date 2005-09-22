// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQED_H
#define HERWIG_ShowerAlphaQED_H
//
// This is the declaration of the ShowerAlphaQED class.

#include "ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value(scale).
 *
 *  @see ShowerAlpha
 */
class ShowerAlphaQED: public ShowerAlpha {

public:

  /**
   * Standard ctors and dtor.
   */
  inline ShowerAlphaQED();
  inline ShowerAlphaQED(const ShowerAlphaQED &);
  virtual ~ShowerAlphaQED();

  /**
   * The following two methods are equivalent to the QCD ones 
   * and are necessary to make use of the virtuality of ShowerAlpha
   * at other places. 
   */

  /**
   * It returns the running alpha value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double value(const Energy2 scale);

  /**
   * It returns the running alpha value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  inline virtual double overestimateValue();

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

private:

  /**
   * Describe an concrete class with persistent data.
   */
  static ClassDescription<ShowerAlphaQED> initShowerAlphaQED;

  /**
   * Private and non-existent assignment operator.
   */
  ShowerAlphaQED & operator=(const ShowerAlphaQED &);

private: 

  double _alpha; 
  
};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerAlphaQED.
 */
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQED,1> {
  typedef Herwig::ShowerAlpha NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ShowerAlphaQED>: 
    public ClassTraitsBase<Herwig::ShowerAlphaQED> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ShowerAlphaQED"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "ShowerAlphaQED.icc"

#endif /* HERWIG_ShowerAlphaQCD_H */
