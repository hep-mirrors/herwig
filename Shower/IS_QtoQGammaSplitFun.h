// -*- C++ -*-
#ifndef HERWIG_IS_QtoQGammaSplitFun_H
#define HERWIG_IS_QtoQGammaSplitFun_H
//
// This is the declaration of the <!id>IS_QtoQGammaSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the exact Next-to-Leading-Order (NLO) <BR>
// Initial State splitting function <I>Q-&GT;QGamma</I>. <BR>
// If you want to use instead the Leading-Order (LO) one, then <BR> 
// do *not* define the virtual methods below.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:QtoQGammaSplitFun.html">QtoQGammaSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGammaSplitFun.html">FS_QtoQGammaSplitFun.h</a>.
// 

#include "QtoQGammaSplitFun.h"


namespace Herwig {

using namespace ThePEG;

class IS_QtoQGammaSplitFun: public QtoQGammaSplitFun {

public:

  inline IS_QtoQGammaSplitFun();
  inline IS_QtoQGammaSplitFun(const IS_QtoQGammaSplitFun &);
  virtual ~IS_QtoQGammaSplitFun();
  // Standard ctors and dtor.

  inline IS_QtoQGammaSplitFun( const long inputIdQuark );
  inline IS_QtoQGammaSplitFun( const long inputIdQuark, const Energy inputMassQuark );

  // virtual Complex fullFun( const double z, const double phi );
  // virtual Complex integratedFun( const double z );
  // virtual Complex fullFunWithHelicities( const double z, const double phi,
  //    		 		    const int h0, const int h1, const int h2 );
  // virtual Complex integratedFunWithHelicities( const double z,
  //	        				  const int h0, const int h1, const int h2 );
  // These methods should be defined only if you want to
  // use the exact Next-to-Leading-Order (NLO) values of the 
  // Initial State splitting function <I>Q-&GT;QGamma</I>, evaluated in terms of 
  // some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and
  // helicities of the three particles.
  // Notice that if you are happy with the LO splitting function, then
  // you should *not* override the virtual methods defined in the
  // base class <!class>QtoQGammaSplitFun<!!class>.

private:

  IS_QtoQGammaSplitFun & operator=(const IS_QtoQGammaSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_QtoQGammaSplitFun.icc"

#endif /* HERWIG_IS_QtoQGammaSplitFun_H */
