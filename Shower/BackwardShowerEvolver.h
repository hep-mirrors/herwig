// -*- C++ -*-
#ifndef HERWIG_BackwardShowerEvolver_H
#define HERWIG_BackwardShowerEvolver_H
//
// This is the declaration of the <!id>BackwardShowerEvolver<!!id> class.
//
// This class is responsible for the backward evolution of a space-like particles <BR>
// (and recursively to all their time-like radiation products). <BR> 
//
// ***LOOKHERE*** At the moment at least, the class has no member data
//                and it has only one method. So it could be easily
//                replaced with a method in <!id>InsideRangeShowerEvolver<!!id>.
//                The same of course for <!id>ForwardShowerEvolver<!!id> class.
//                However, we could end up with a too giant class
//                <!id>InsideRangeShowerEvolver<!!id>, so I prefer to keep
//                the classes <!id>ForwardShowerEvolver<!!id> and <!id>BackwardShowerEvolver<!!id>.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>, <BR>
// <a href="http:RhoDMatrixPropagator.html">RhoDMatrixPropagator.h</a>, <BR>
// <a href="http:ForwardShowerEvolver.html">ForwardShowerEvolver.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "Pythia7/Handlers/PartialCollisionHandler.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "ShowerConfig.h"
#include "SplittingGenerator.h"
#include "RhoDMatrixPropagator.h"
#include "ForwardShowerEvolver.h"


namespace Herwig {

using namespace Pythia7;

class BackwardShowerEvolver: public Pythia7::HandlerBase {

public:

  inline BackwardShowerEvolver();
  inline BackwardShowerEvolver(const BackwardShowerEvolver &);
  virtual ~BackwardShowerEvolver();
  // Standard ctors and dtor.

  bool spaceLikeShower( tPartCollHdlPtr ch, 
		        const tShoConstrPtr showerConstrainer, 
		        const tMECorrectionPtr meCorrectionPtr,
		        tShoParPtr particle, CollecShoParPtr & collecShoPar ) 
    throw (Veto, Stop, Exception);
  // It does the backward evolution of the space-like input <!id>particle<!!id> 
  // (and recursively for all its time-like radiation products).
  // accepting only emissions which conforms to the <!id>showerConstrainer<!!id>
  // and soft matrix element correction pointed by <!id>meCorrectionPtr<!!id>.
  // The <!id>ParticleCollisionHandler<!!id> object is needed to access the PDFs.
  // If at least one emission has occurred then the method returns true
  // and all the new created <!id>ShowerParticle<!!id> objects (but not the input
  // particle) are added to the collection <!id>collecShoPar<!!id> (which can
  // contain, at the beginning of the method, either the full collection
  // of <!id>ShowerParticle<!!id> already created so far by the showering, 
  // or being empty: the choice is up to the caller).  

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

  static ClassDescription<BackwardShowerEvolver> initBackwardShowerEvolver;
  // Describe a concrete class with persistent data.

  BackwardShowerEvolver & operator=(const BackwardShowerEvolver &);
  //  Private and non-existent assignment operator.

  Ptr<SplittingGenerator>::pointer _pointerSplittingGenerator;
  Ptr<RhoDMatrixPropagator>::pointer _pointerRhoDMatrixPropagator;
  Ptr<ForwardShowerEvolver>::pointer _pointerForwardShowerEvolver;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of BackwardShowerEvolver.
template <>
struct BaseClassTrait<Herwig::BackwardShowerEvolver,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::BackwardShowerEvolver>: public ClassTraitsBase<Herwig::BackwardShowerEvolver> {
  static string className() { return "/Herwig++/BackwardShowerEvolver"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "BackwardShowerEvolver.icc"

#endif /* HERWIG_BackwardShowerEvolver_H */
