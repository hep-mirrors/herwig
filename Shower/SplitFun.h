// -*- C++ -*-
#ifndef HERWIG_SplitFun_H
#define HERWIG_SplitFun_H
//
// This is the declaration of the <!id>SplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is abstract class from which all splitting function classes, <BR>
// whatever their interaction type (QCD, QED, EWK,...) and multiplicity <BR>
// <I>1-&GT;N</I>, are inheriting from.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:SplitFun1to3.html">SplitFun1to3.h</a>.

#include "Pythia7/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "ShowerIndex.h"


namespace Herwig {

using namespace Pythia7;

class SplitFun: public Pythia7::HandlerBase {

public:

  inline SplitFun();
  inline SplitFun(const SplitFun &);
  virtual ~SplitFun();
  // Standard ctors and dtor.

  inline SplitFun( const ShowerIndex::InteractionType interaction,
		   const int inputNumBranchingProducts,
                   const long inputIdEmitter, const Energy inputMassEmitter);
  // Specifies the interaction type of the vertex <I>0-&GT;1+2+...+N</I> 
  // the number of branching products (<I>N</I>), and the PDG id and mass of the emitter.
  // (Notice that the id and masses of branching products are of 
  // responsability of the classes that inherit from this class).

  inline ShowerIndex::InteractionType interactionType() const;
  // Type (QCD, QED, EWK,...) of the emission (branching) vertex.
  
  inline int numBranchingProducts() const;
  // The number of branching products.

  inline long idEmitter() const;
  inline Energy massEmitter() const;
  // PDG id and mass of the emitting (showering) particle.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static AbstractClassDescription<SplitFun> initSplitFun;
  // Describe an abstract base class with persistent data.

  SplitFun & operator=(const SplitFun &);
  //  Private and non-existent assignment operator.
  
  ShowerIndex::InteractionType _interaction;
  int _numProducts;
  long _idEmitter;
  Energy _mEmitter;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of SplitFun.
template <>
struct BaseClassTrait<Herwig::SplitFun,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::SplitFun>: public ClassTraitsBase<Herwig::SplitFun> {
  static string className() { return "/Herwig++/SplitFun"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "SplitFun.icc"

#endif /* HERWIG_SplitFun_H */
