// -*- C++ -*-
#ifndef HERWIG_IS_GtoGGSplitFun_H
#define HERWIG_IS_GtoGGSplitFun_H
//
// This is the declaration of the <!id>IS_GtoGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the exact Next-to-Leading-Order (NLO) <BR>
// Initial State splitting function <I>G-&GT;GG</I>. <BR> 
// If you want to use instead the Leading-Order (LO) one, then <BR> 
// do *not* define the virtual methods below.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GtoGGSplitFun.html">GtoGGSplitFun.h</a>, <BR>
// <a href="http:FS_GtoGGSplitFun.html">FS_GtoGGSplitFun.h</a>.
// 

#include "GtoGGSplitFun.h"


namespace Herwig {

using namespace Pythia7;

class IS_GtoGGSplitFun: public GtoGGSplitFun {

public:

  inline IS_GtoGGSplitFun();
  inline IS_GtoGGSplitFun(const IS_GtoGGSplitFun &);
  virtual ~IS_GtoGGSplitFun();
  // Standard ctors and dtor.

  // virtual Complex fullFun( const double z, const double phi );
  // virtual Complex integratedFun( const double z );
  // virtual Complex fullFunWithHelicities( const double z, const double phi,
  //	   				    const int h0, const int h1, const int h2 );
  // virtual Complex integratedFunWithHelicities( const double z,
  //					          const int h0, const int h1, const int h2 );
  // These methods should be defined only if you want to
  // use the exact Next-to-Leading-Order (NLO) values of the 
  // Initial State splitting function <I>G-&GT;GG</I>, evaluated in terms of 
  // some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.
  // Notice that if you are happy with the LO splitting function, then
  // you should *not* override the virtual methods defined in the
  // base class <!class>GtoGGSplitFun<!!class>.

private:

  IS_GtoGGSplitFun & operator=(const IS_GtoGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_GtoGGSplitFun.icc"

#endif /* HERWIG_IS_GtoGGSplitFun_H */
