// -*- C++ -*-
#ifndef HERWIG_IS_GtoQQbarSplitFun_H
#define HERWIG_IS_GtoQQbarSplitFun_H
//
// This is the declaration of the <!id>IS_GtoQQbarSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the exact Next-to-Leading-Order (NLO) <BR>
// Initial State splitting function <I>G-&GT;QQbar</I>. <BR> 
// If you want to use instead the Leading-Order (LO) one, then <BR> 
// do *not* define the virtual methods below.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GtoQQbarSplitFun.html">GtoQQbarSplitFun.h</a>, <BR>
// <a href="http:FS_GtoQQbarSplitFun.html">FS_GtoQQbarSplitFun.h</a>.
// 

#include "GtoQQbarSplitFun.h"


namespace Herwig {

using namespace ThePEG;

class IS_GtoQQbarSplitFun: public GtoQQbarSplitFun {

public:

  inline IS_GtoQQbarSplitFun();
  inline IS_GtoQQbarSplitFun(const IS_GtoQQbarSplitFun &);
  virtual ~IS_GtoQQbarSplitFun();
  // Standard ctors and dtor.

  inline IS_GtoQQbarSplitFun( const long inputIdQuark );
  inline IS_GtoQQbarSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  // virtual Complex fullFun( const double z, const double phi );
  // virtual Complex integratedFun( const double z );
  // virtual Complex fullFunWithHelicities( const double z, const double phi,
  //			         	    const int h0, const int h1, const int h2 );
  // virtual Complex integratedFunWithHelicities( const double z,
  //					          const int h0, const int h1, const int h2 );
  // These methods should be defined only if you want to
  // use the exact Next-to-Leading-Order (NLO) values of the 
  // Initial State splitting function <I>G-&GT;QQbar</I>, evaluated in terms of 
  // some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and
  // helicities of the three particles.
  // Notice that if you are happy with the LO splitting function, then
  // you should *not* override the virtual methods defined in the
  // base class <!class>GtoQQbarSplitFun<!!class>.

private:

  IS_GtoQQbarSplitFun & operator=(const IS_GtoQQbarSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_GtoQQbarSplitFun.icc"

#endif /* HERWIG_IS_GtoQQbarSplitFun_H */
