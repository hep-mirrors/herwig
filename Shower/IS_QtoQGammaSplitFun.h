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

using namespace Pythia7;

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

  static ClassDescription<IS_QtoQGammaSplitFun> initIS_QtoQGammaSplitFun;
  // Describe an abstract base class with persistent data.

  IS_QtoQGammaSplitFun & operator=(const IS_QtoQGammaSplitFun &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of IS_QtoQGammaSplitFun.
template <>
struct BaseClassTrait<Herwig::IS_QtoQGammaSplitFun,1> {
  typedef Herwig::QtoQGammaSplitFun NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::IS_QtoQGammaSplitFun>: public ClassTraitsBase<Herwig::IS_QtoQGammaSplitFun> {
  static string className() { return "/Herwig++/IS_QtoQGammaSplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "IS_QtoQGammaSplitFun.icc"

#endif /* HERWIG_IS_QtoQGammaSplitFun_H */
