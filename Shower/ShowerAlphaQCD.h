// -*- C++ -*-
#ifndef HERWIG_ShowerAlphaQCD_H
#define HERWIG_ShowerAlphaQCD_H
//
// This is the declaration of the <!id>ShowerAlphaQCD<!!id> class.
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
#include "ShowerIndex.h"
namespace Herwig {

using namespace Pythia7;

class ShowerAlphaQCD: public ShowerAlpha {

public:

  inline ShowerAlphaQCD();
  inline ShowerAlphaQCD(const ShowerAlphaQCD &);
  virtual ~ShowerAlphaQCD();
  // Standard ctors and dtor.

  virtual double value(const Energy2 scale);
  // It returns the running alpha value evaluated at the input scale
  // multiplied by the scale factor <!id>scaleFactor()<!!id>.

  virtual double overestimateValue();
  // It returns the running alpha value evaluated at the input scale
  // multiplied by the scale factor <!id>scaleFactor()<!!id>.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

  void persistentOutput(PersistentOStream & os) const;
  void persistentInput(PersistentIStream & is, int);

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

private:

  static ClassDescription<ShowerAlphaQCD> initShowerAlphaQCD;
  // Describe an concrete class with persistent data.

  ShowerAlphaQCD & operator=(const ShowerAlphaQCD &);
  //  Private and non-existent assignment operator.

private: 

  double alphaTwoLoop(Energy q, Energy lam, short nf); 
  // the two-loop parametrization of alpha_s

  pair<short, Energy> getLamNfTwoLoop(Energy q); 
  // hacked in masses by hand for the moment before proper
  // interfacing...  obtained lambda solutions numerically in
  // Mathematica with my alphas.m using two-loop alphas from PDT2002
  // and as(M_Z=91.2GeV) = 0.118 *** ACHTUNG! *** this HAS to be done
  // automatically acc to the masses and as(M_Z) given by the PDT
  // class (which is supposed to be up-to-date)

  double alpha_s(double q2, double q2min, int type); 
  // A toy parametrization of alpha_s with different parametrizations
  // of the IR behaviour, below <!id>q2min<!!id>, set by
  // <!id>type<!!id>.  Default is <!id>type = 1<!!id>,
  // i.e. <!id>alpha_s<!!id> = 0 below <!id>q2min<!!id>.

  int _asType;
  Energy _Qmin; 

  void test(int i);
};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ShowerAlphaQCD.
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQCD,1> {
  typedef Herwig::ShowerAlpha NthBase;
};



// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ShowerAlphaQCD>: public ClassTraitsBase<Herwig::ShowerAlphaQCD> {
  static string className() { return "/Herwig++/ShowerAlphaQCD"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ShowerAlphaQCD.icc"

#endif /* HERWIG_ShowerAlphaQCD_H */
