// -*- C++ -*-
#ifndef HERWIG_FS_GtoGGGSplitFun_H
#define HERWIG_FS_GtoGGGSplitFun_H
//
// This is the declaration of the <!id>FS_GtoGGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the (concrete) class which describes the <BR>
// Final State <I>G -&GT; G+G+G</I> splitting function.
//
// ***LOOKHERE***  This class is currently kept empty; however, 
//                 if you have to implement it, you should proceed 
//                 similarly to the <!class>FS_GtoGGSplitFun<!!class> class.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GtoGGGSplitFun.html">GtoGGGSplitFun.h</a>, <BR>
// <a href="http:IS_GtoGGGSplitFun.html">IS_GtoGGGSplitFun.h</a>.
// 

#include "GtoGGGSplitFun.h"


namespace Herwig {

using namespace Pythia7;

class FS_GtoGGGSplitFun: public GtoGGGSplitFun {

public:

  inline FS_GtoGGGSplitFun();
  inline FS_GtoGGGSplitFun(const FS_GtoGGGSplitFun &);
  virtual ~FS_GtoGGGSplitFun();
  // Standard ctors and dtor.

private:

  FS_GtoGGGSplitFun & operator=(const FS_GtoGGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "FS_GtoGGGSplitFun.icc"

#endif /* HERWIG_FS_GtoGGGSplitFun_H */
