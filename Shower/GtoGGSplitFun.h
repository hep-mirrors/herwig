// -*- C++ -*-
#ifndef HERWIG_GtoGGSplitFun_H
#define HERWIG_GtoGGSplitFun_H
//
// This is the declaration of the <!id>GtoGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting <BR>
// function for <I>G-&GT;GG</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_GtoGGSplitFun.html">IS_GtoGGSplitFun.h</a>, <BR>
// <a href="http:FS_GtoGGSplitFun.html">FS_GtoGGSplitFun.h</a>.
//  

#include "SplitFun1to2.h"


namespace Herwig {

using namespace ThePEG;

class GtoGGSplitFun: public SplitFun1to2 {

public:

  inline GtoGGSplitFun();
  inline GtoGGSplitFun(const GtoGGSplitFun &);
  virtual ~GtoGGSplitFun();
  // Standard ctors and dtor.

  virtual double fullFun( const double z, const Energy2 qtilde2, const double phi );
  virtual double integratedFun( const double z, const Energy2 qtilde2);
  virtual double fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 );
  virtual double integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 );
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>G-&GT;GG</I>
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

private:

  GtoGGSplitFun & operator=(const GtoGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "GtoGGSplitFun.icc"

#endif /* HERWIG_GtoGGSplitFun_H */
