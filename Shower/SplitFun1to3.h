// -*- C++ -*-
#ifndef HERWIG_SplitFun1to3_H
#define HERWIG_SplitFun1to3_H
//
// This is the declaration of the <!id>SplitFun1to3<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is an abstract class which describes the common interface <BR>
// of all <I>1-&GT;3</I> splitting functions, both for Initial State <BR>
// and Final State Radiation. 
//
// ***LOOKHERE***  This class is currently kept empty; however, 
//                 if you have to implement it, you should proceed 
//                 similarly to the <!class>SplitFun1to2<!!class> class.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun.html">SplitFun.h</a>, <BR>
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:GtoGGGSplitFun.html">GtoGGGSplitFun.h</a>.
// 

#include "SplitFun.h"


namespace Herwig {

using namespace Pythia7;

class SplitFun1to3: public SplitFun {

public:

  inline SplitFun1to3();
  inline SplitFun1to3(const SplitFun1to3 &);
  virtual ~SplitFun1to3();
  // Standard ctors and dtor.

private:

  SplitFun1to3 & operator=(const SplitFun1to3 &);
  //  Private and non-existent assignment operator.

};

}

#include "SplitFun1to3.icc"

#endif /* HERWIG_SplitFun1to3_H */
