// -*- C++ -*-
#ifndef HERWIG_IS_QtoQGSplitFun_H
#define HERWIG_IS_QtoQGSplitFun_H
//
// This is the declaration of the <!id>IS_QtoQGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the exact Next-to-Leading-Order (NLO) <BR>
// Initial State splitting function <I>Q-&GT;QG</I>. <BR> 
// If you want to use instead the Leading-Order (LO) one, then <BR> 
// do *not* define the virtual methods below.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:QtoQGSplitFun.html">QtoQGSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGSplitFun.html">FS_QtoQGSplitFun.h</a>.
// 

#include "QtoQGSplitFun.h"


namespace Herwig {

using namespace ThePEG;

class IS_QtoQGSplitFun: public QtoQGSplitFun {

public:

  inline IS_QtoQGSplitFun();
  inline IS_QtoQGSplitFun(const IS_QtoQGSplitFun &);
  virtual ~IS_QtoQGSplitFun();
  // Standard ctors and dtor.

  inline IS_QtoQGSplitFun( const long inputIdQuark );
  inline IS_QtoQGSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  // virtual Complex fullFun( const double z, const double phi );
  // virtual Complex integratedFun( const double z );
  // virtual Complex fullFunWithHelicities( const double z, const double phi,
  //	      				    const int h0, const int h1, const int h2 );
  // virtual Complex integratedFunWithHelicities( const double z,
  //					          const int h0, const int h1, const int h2 );
  // These methods should be defined only if you want to
  // use the exact Next-to-Leading-Order (NLO) values of the 
  // Initial State splitting function <I>Q-&GT;QG</I>, evaluated in terms of 
  // some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.
  // Notice that if you are happy with the LO splitting function, then
  // you should *not* override the virtual methods defined in the
  // base class <!class>QtoQGSplitFun<!!class>.

private:

  IS_QtoQGSplitFun & operator=(const IS_QtoQGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_QtoQGSplitFun.icc"

#endif /* HERWIG_IS_QtoQGSplitFun_H */
