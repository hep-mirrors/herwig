// -*- C++ -*-
#ifndef HERWIG_GtoGGGSplitFun_H
#define HERWIG_GtoGGGSplitFun_H
//
// This is the declaration of the <!id>GtoGGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is an abstract class which describes the common <BR>
// part of the <I>G-&GT;GGG</I> splitting function for both Initial State <BR> 
// and Final State Radiation. It should factorize the common <BR>
// expression at LO, a part the coupling constant which is <BR>
// kept separated for I.S.R. and F.S.R. for systematics evaluation.
//
// ***LOOKHERE***  This class is currently kept empty; however, 
//                 if you have to implement it, you should proceed 
//                 similarly to the <!class>GtoGGSplitFun<!!class> class.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to3.html">SplitFun1to3.h</a>, <BR>
// <a href="http:IS_GtoGGGSplitFun.html">IS_GtoGGGSplitFun.h</a>, <BR>
// <a href="http:FS_GtoGGGSplitFun.html">FS_GtoGGGSplitFun.h</a>.
// 

#include "SplitFun1to3.h"


namespace Herwig {

using namespace ThePEG;

class GtoGGGSplitFun: public SplitFun1to3 {

public:

  inline GtoGGGSplitFun();
  inline GtoGGGSplitFun(const GtoGGGSplitFun &);
  virtual ~GtoGGGSplitFun();
  // Standard ctors and dtor.

private:

  GtoGGGSplitFun & operator=(const GtoGGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

#include "GtoGGGSplitFun.icc"

#endif /* HERWIG_GtoGGGSplitFun_H */
