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

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

private:

  static ClassDescription<IS_GtoGGGSplitFun> initIS_GtoGGGSplitFun;
  // Describe an abstract base class with persistent data.

  IS_GtoGGGSplitFun & operator=(const IS_GtoGGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of IS_GtoGGGSplitFun.
template <>
struct BaseClassTrait<Herwig::IS_GtoGGGSplitFun,1> {
  typedef Herwig::GtoGGGSplitFun NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::IS_GtoGGGSplitFun>: public ClassTraitsBase<Herwig::IS_GtoGGGSplitFun> {
  static string className() { return "/Herwig++/IS_GtoGGGSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "IS_GtoGGGSplitFun.icc"

#endif /* HERWIG_IS_GtoGGGSplitFun_H */
