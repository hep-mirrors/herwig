// -*- C++ -*-
#ifndef HERWIG_QtoQGammaSplitFun_H
#define HERWIG_QtoQGammaSplitFun_H
//
// This is the declaration of the <!id>QtoQGammaSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting
// function for <I>Q-&GT;QGamma</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_QtoQGammaSplitFun.html">IS_QtoQGammaSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGammaSplitFun.html">FS_QtoQGammaSplitFun.h</a>.
//  

#include "SplitFun1to2.h"


namespace Herwig {

using namespace Pythia7;

class QtoQGammaSplitFun: public SplitFun1to2 {

public:

  inline QtoQGammaSplitFun();
  inline QtoQGammaSplitFun(const QtoQGammaSplitFun &);
  virtual ~QtoQGammaSplitFun();
  // Standard ctors and dtor.

  inline QtoQGammaSplitFun( const long inputIdQuark );
  inline QtoQGammaSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  virtual Complex fullFun( const double z, const Energy2 qtilde2, const double phi );
  virtual Complex integratedFun( const double z, const Energy2 qtilde2);
  virtual Complex fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 );
  virtual Complex integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 );
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>Q-&GT;QGamma</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.

  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );
  // See long comment on this method on class <!class>SplitFun1to2<!!class>.
  // Remember that the first branching product is considered the quark
  // and the second one the photon.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static AbstractClassDescription<QtoQGammaSplitFun> initQtoQGammaSplitFun;
  // Describe an abstract base class with persistent data.

  QtoQGammaSplitFun & operator=(const QtoQGammaSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of QtoQGammaSplitFun.
template <>
struct BaseClassTrait<Herwig::QtoQGammaSplitFun,1> {
  typedef Herwig::SplitFun1to2 NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::QtoQGammaSplitFun>: public ClassTraitsBase<Herwig::QtoQGammaSplitFun> {
  static string className() { return "/Herwig++/QtoQGammaSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "QtoQGammaSplitFun.icc"

#endif /* HERWIG_QtoQGammaSplitFun_H */
