// -*- C++ -*-
#ifndef HERWIG_IS_GtoGGGSplitFun_H
#define HERWIG_IS_GtoGGGSplitFun_H
//
// This is the declaration of the <!id>IS_GtoGGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the (concrete) class which describes the <BR>
// Initial State <I>G -&GT; G+G+G</I> splitting function.
//
// ***LOOKHERE***  This class is currently kept empty; however, 
//                 if you have to implement it, you should proceed 
//                 similarly to the <!class>IS_GtoGGSplitFun<!!class> class.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GtoGGGSplitFun.html">GtoGGGSplitFun.h</a>, <BR>
// <a href="http:FS_GtoGGGSplitFun.html">FS_GtoGGGSplitFun.h</a>.
// 

#include "GtoGGGSplitFun.h"


namespace Herwig {

using namespace Pythia7;

class IS_GtoGGGSplitFun: public GtoGGGSplitFun {

public:

  inline IS_GtoGGGSplitFun();
  inline IS_GtoGGGSplitFun(const IS_GtoGGGSplitFun &);
  virtual ~IS_GtoGGGSplitFun();
  // Standard ctors and dtor.

private:

  IS_GtoGGGSplitFun & operator=(const IS_GtoGGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_GtoGGGSplitFun.icc"

#endif /* HERWIG_IS_GtoGGGSplitFun_H */
