// -*- C++ -*-
#ifndef HERWIG_QtoQGSplitFun_H
#define HERWIG_QtoQGSplitFun_H
//
// This is the declaration of the <!id>QtoQGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting <BR>
// function for <I>Q-&GT;QG</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_QtoQGSplitFun.html">IS_QtoQGSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGSplitFun.html">FS_QtoQGSplitFun.h</a>.
// 

#include "SplitFun1to2.h"


namespace Herwig {

using namespace Pythia7;

class QtoQGSplitFun: public SplitFun1to2 {

public:

  inline QtoQGSplitFun();
  inline QtoQGSplitFun(const QtoQGSplitFun &);
  virtual ~QtoQGSplitFun();
  // Standard ctors and dtor.

  inline QtoQGSplitFun( const long inputIdQuark );
  inline QtoQGSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  virtual Complex fullFun( const double z, const Energy2 qtilde2, const double phi );
  virtual Complex integratedFun( const double z, const Energy2 qtilde2 );
  virtual Complex fullFunWithHelicities( const double z, const Energy2 qtilde2,
					 const double phi, 
					 const int h0, const int h1, const int h2 );
  virtual Complex integratedFunWithHelicities( const double z, 
					       const Energy2 qtilde2, 
					       const int h0, const int h1, 
					       const int h2 );
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>Q-&GT;QG</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.

  virtual Complex overestimateIntegratedFun( const double z );
  virtual Complex integOverIntegratedFun(const double z); 
  virtual Complex invIntegOverIntegratedFun(const double r); 

  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );
  // See long comment on this method on class <!class>SplitFun1to2<!!class>.
  // Remember that the first branching product is considered the quark
  // and the second one the gluon.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static AbstractClassDescription<QtoQGSplitFun> initQtoQGSplitFun;
  // Describe an abstract base class with persistent data.

  QtoQGSplitFun & operator=(const QtoQGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of QtoQGSplitFun.
template <>
struct BaseClassTrait<Herwig::QtoQGSplitFun,1> {
  typedef Herwig::SplitFun1to2 NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::QtoQGSplitFun>: public ClassTraitsBase<Herwig::QtoQGSplitFun> {
  static string className() { return "/Herwig++/QtoQGSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "QtoQGSplitFun.icc"

#endif /* HERWIG_QtoQGSplitFun_H */
