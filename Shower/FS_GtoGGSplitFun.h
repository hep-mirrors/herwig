// -*- C++ -*-
#ifndef HERWIG_FS_GtoGGSplitFun_H
#define HERWIG_FS_GtoGGSplitFun_H
//
// This is the declaration of the <!id>FS_GtoGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the exact Next-to-Leading-Order (NLO) <BR>
// Final State splitting function <I>G-&GT;GG</I>. <BR> 
// If you want to use instead the Leading-Order (LO) one, then <BR> 
// do *not* define the virtual methods below.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GtoGGSplitFun.html">GtoGGSplitFun.h</a>, <BR>
// <a href="http:IS_GtoGGSplitFun.html">IS_GtoGGSplitFun.h</a>.
// 

#include "GtoGGSplitFun.h"


namespace Herwig {

using namespace Pythia7;

class FS_GtoGGSplitFun: public GtoGGSplitFun {

public:

  inline FS_GtoGGSplitFun();
  inline FS_GtoGGSplitFun(const FS_GtoGGSplitFun &);
  virtual ~FS_GtoGGSplitFun();
  // Standard ctors and dtor.

  // virtual Complex fullFun( const double z, const double phi );
  // virtual Complex integratedFun( const double z );
  // virtual Complex fullFunWithHelicities( const double z, const double phi,
  //					    const int h0, const int h1, const int h2 );
  // virtual Complex integratedFunWithHelicities( const double z,
  //					          const int h0, const int h1, const int h2 );
  // These methods should be defined only if you want to
  // use the exact Next-to-Leading-Order (NLO) values of the 
  // Final State splitting function <I>G-&GT;GG</I>, evaluated in terms of 
  // some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.
  // Notice that if you are happy with the LO splitting function, then
  // you should *not* override the virtual methods defined in the
  // base class <!class>GtoGGSplitFun<!!class>.

private:

  FS_GtoGGSplitFun & operator=(const FS_GtoGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "FS_GtoGGSplitFun.icc"

#endif /* HERWIG_FS_GtoGGSplitFun_H */
