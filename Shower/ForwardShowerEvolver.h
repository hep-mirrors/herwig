// -*- C++ -*-
#ifndef HERWIG_ForwardShowerEvolver_H
#define HERWIG_ForwardShowerEvolver_H
//
// This is the declaration of the <!id>ForwardShowerEvolver<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for the forward evolution of a time-like particles <BR>
// (and recursively to all their radiation products). <BR>
// It also treats the special case of the showering of a time-like decaying particle, <BR>
// in which the emissions have reversed angular ordering. <BR>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>, <BR>
// <a href="http:RhoDMatrixPropagator.html">RhoDMatrixPropagator.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "SplittingGenerator.h"
#include "RhoDMatrixPropagator.h"


namespace Herwig {

using namespace ThePEG;

class ForwardShowerEvolver: public ThePEG::HandlerBase {

public:

  inline ForwardShowerEvolver();
  inline ForwardShowerEvolver(const ForwardShowerEvolver &);
  virtual ~ForwardShowerEvolver();
  // Standard ctors and dtor.

  bool timeLikeShower( tPartCollHdlPtr ch, 
		       const tShoConstrPtr showerConstrainer, 
		       const tMECorrectionPtr meCorrectionPtr,
		       tShowerParticlePtr particle, 
		       ShowerParticleVector & collecShoPar,
		       const bool specialDecay = false ) throw (Veto, Stop, Exception); 
  // It does the forward evolution of the time-like input <!id>particle<!!id>
  // (and recursively for all its radiation products).
  // accepting only emissions which conforms to the <!id>showerConstrainer<!!id>
  // and soft matrix element correction pointed by  <!id>meCorrectionPtr<!!id>.
  // In the case that <!id>specialDecay<!!id> is true then the forward evolution
  // is done with reverse angular ordering, as it should be for radiation
  // emitted by a decaying particle. 
  // If at least one emission has occurred then the method returns true
  // and all the new created <!id>ShowerParticle<!!id> objects (but not the input 
  // particle) are added to the collection <!id>collecShoPar<!!id> (which can
  // contain, at the beginning of the method, either the full collection
  // of <!id>ShowerParticle<!!id> already created so far by the showering, 
  // or being empty: the choice is up to the caller).  

private:
  bool MEVeto(tcPPtr, const Energy &, const double &);

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

  static ClassDescription<ForwardShowerEvolver> initForwardShowerEvolver;
  // Describe a concrete class with persistent data.

  ForwardShowerEvolver & operator=(const ForwardShowerEvolver &);
  //  Private and non-existent assignment operator.

  Ptr<SplittingGenerator>::pointer _pointerSplittingGenerator;
  Ptr<RhoDMatrixPropagator>::pointer _pointerRhoDMatrixPropagator;

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ForwardShowerEvolver.
template <>
struct BaseClassTrait<Herwig::ForwardShowerEvolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ForwardShowerEvolver>: public ClassTraitsBase<Herwig::ForwardShowerEvolver> {
  static string className() { return "/Herwig++/ForwardShowerEvolver"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ForwardShowerEvolver.icc"

#endif /* HERWIG_ForwardShowerEvolver_H */
