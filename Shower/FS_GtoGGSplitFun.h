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

  static ClassDescription<FS_GtoGGSplitFun> initFS_GtoGGSplitFun;
  // Describe an abstract base class with persistent data.

  FS_GtoGGSplitFun & operator=(const FS_GtoGGSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of FS_GtoGGSplitFun.
template <>
struct BaseClassTrait<Herwig::FS_GtoGGSplitFun,1> {
  typedef Herwig::GtoGGSplitFun NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::FS_GtoGGSplitFun>: public ClassTraitsBase<Herwig::FS_GtoGGSplitFun> {
  static string className() { return "/Herwig++/FS_GtoGGSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "FS_GtoGGSplitFun.icc"

#endif /* HERWIG_FS_GtoGGSplitFun_H */
