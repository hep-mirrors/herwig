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

  inline QtoQGSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  virtual double fullFun( const double z, const Energy2 qtilde2, const double phi );
  virtual double integratedFun( const double z, const Energy2 qtilde2 );
  virtual double fullFunWithHelicities( const double z, const Energy2 qtilde2,
					 const double phi, 
					 const int h0, const int h1, const int h2 );
  virtual double integratedFunWithHelicities( const double z, 
					       const Energy2 qtilde2, 
					       const int h0, const int h1, 
					       const int h2 );
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>Q-&GT;QG</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.

  virtual double overestimateIntegratedFun( const double z );
  virtual double integOverIntegratedFun(const double z); 
  virtual double invIntegOverIntegratedFun(const double r); 

  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );
  // See long comment on this method on class <!class>SplitFun1to2<!!class>.
  // Remember that the first branching product is considered the quark
  // and the second one the gluon.

private:

  QtoQGSplitFun & operator=(const QtoQGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "QtoQGSplitFun.icc"

#endif /* HERWIG_QtoQGSplitFun_H */
