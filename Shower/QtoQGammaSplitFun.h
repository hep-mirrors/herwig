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

private:

  QtoQGammaSplitFun & operator=(const QtoQGammaSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "QtoQGammaSplitFun.icc"

#endif /* HERWIG_QtoQGammaSplitFun_H */
