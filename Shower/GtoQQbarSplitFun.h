// -*- C++ -*-
#ifndef HERWIG_GtoQQbarSplitFun_H
#define HERWIG_GtoQQbarSplitFun_H
//
// This is the declaration of the <!id>GtoQQbarSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting
// function for <I>G-&GT;QQbar</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_GtoQQbarSplitFun.html">IS_GtoQQbarSplitFun.h</a>, <BR>
// <a href="http:FS_GtoQQbarSplitFun.html">FS_GtoQQbarSplitFun.h</a>.
//  

#include "SplitFun1to2.h"


namespace Herwig {

using namespace Pythia7;

class GtoQQbarSplitFun: public SplitFun1to2 {

public:

  inline GtoQQbarSplitFun();
  inline GtoQQbarSplitFun(const GtoQQbarSplitFun &);
  virtual ~GtoQQbarSplitFun();
  // Standard ctors and dtor.

  inline GtoQQbarSplitFun( const long inputIdQuark );
  inline GtoQQbarSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  virtual Complex fullFun( const double z, const Energy2 qtilde2, const double phi );
  virtual Complex integratedFun( const double z, const Energy2 qtilde2 );
  virtual Complex fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 );
  virtual Complex integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 );
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>G-&GT;QQbar</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and
  // helicities of the three particles.

  virtual Complex overestimateIntegratedFun( const double z );
  virtual Complex integOverIntegratedFun(const double z); 
  virtual Complex invIntegOverIntegratedFun(const double r); 


public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static AbstractClassDescription<GtoQQbarSplitFun> initGtoQQbarSplitFun;
  // Describe an abstract base class with persistent data.

  GtoQQbarSplitFun & operator=(const GtoQQbarSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of GtoQQbarSplitFun.
template <>
struct BaseClassTrait<Herwig::GtoQQbarSplitFun,1> {
  typedef Herwig::SplitFun1to2 NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::GtoQQbarSplitFun>: public ClassTraitsBase<Herwig::GtoQQbarSplitFun> {
  static string className() { return "/Herwig++/GtoQQbarSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "GtoQQbarSplitFun.icc"

#endif /* HERWIG_GtoQQbarSplitFun_H */
