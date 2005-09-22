// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQCD_H
#define HERWIG_ShowerAlphaQCD_H
//
// This is the declaration of the ShowerAlphaQCD class.

#include "ShowerAlpha.h"
#include "ShowerIndex.h"
namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value(scale).
 *
 *  @see ShowerAlpha
 */
class ShowerAlphaQCD: public ShowerAlpha {

public:

  /**
   * Standard ctors and dtor.
   */
  inline ShowerAlphaQCD();
  inline ShowerAlphaQCD(const ShowerAlphaQCD &);
  virtual ~ShowerAlphaQCD();

  /**
   * It returns the running alpha value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double value(const Energy2 scale);

  /**
   * It returns the running alpha value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double overestimateValue();

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

  void persistentOutput(PersistentOStream & os) const;
  void persistentInput(PersistentIStream & is, int);

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
 static ClassDescription<ShowerAlphaQCD> initShowerAlphaQCD;

  /**
   * Private and non-existent assignment operator.
   */
  ShowerAlphaQCD & operator=(const ShowerAlphaQCD &);

private: 

  /**
   * The two-loop parametrization of alpha_s. 
   */
  double alphaTwoLoop(Energy q, Energy lam, short nf); 

  /**
   * Hacked in masses by hand for the moment before proper
   * interfacing...  obtained lambda solutions numerically in
   * Mathematica with my alphas.m using two-loop alphas from PDT2002
   * and as(M_Z=91.2GeV) = 0.118 *** ACHTUNG! *** this HAS to be done
   * automatically acc to the masses and as(M_Z) given by the PDT
   * class (which is supposed to be up-to-date).
   */
  pair<short, Energy> getLamNfTwoLoop(Energy q); 

  /**
   * A toy parametrization of alpha_s with different parametrizations
   * of the IR behaviour, below <!id>q2min<!!id>, set by type.  
   * Default is type = 1, i.e. alpha_s = 0 below q2min.
   */
  double alpha_s(double q2, double q2min, int type); 

  int _asType;
  Energy _Qmin; 

  void test(int i);
};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerAlphaQCD.
 */
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQCD,1> {
  typedef Herwig::ShowerAlpha NthBase;
};


/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ShowerAlphaQCD>: 
    public ClassTraitsBase<Herwig::ShowerAlphaQCD> {
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ShowerAlphaQCD"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "ShowerAlphaQCD.icc"

#endif /* HERWIG_ShowerAlphaQCD_H */
