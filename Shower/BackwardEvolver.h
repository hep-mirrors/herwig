// -*- C++ -*-
#ifndef HERWIG_BackwardShowerEvolver_H
#define HERWIG_BackwardShowerEvolver_H
//
// This is the declaration of the <!id>BackwardShowerEvolver<!!id> class.
//
// This class is responsible for the backward evolution of a space-like particles <BR>
// (and recursively to all their time-like radiation products). <BR> 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>, <BR>
// <a href="http:ForwardShowerEvolver.html">ForwardShowerEvolver.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/PartialCollisionHandler.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerConfig.h"
#include "SplittingGenerator.h"
#include "ForwardEvolver.h"


namespace Herwig {

using namespace ThePEG;

class BackwardEvolver: public ThePEG::HandlerBase {

public:

  inline BackwardEvolver();
  inline BackwardEvolver(const BackwardEvolver &);
  virtual ~BackwardEvolver();
  // Standard ctors and dtor.

  bool spaceLikeShower( tPartCollHdlPtr ch, 
		        const tShowerVarsPtr showerVariables, 
		        //const tMECorrectionPtr meCorrectionPtr,
		        tShowerParticlePtr particle, 
			ShowerParticleVector &allShowerParticles) 
    throw (Veto, Stop, Exception);
  // It does the backward evolution of the space-like input <!id>particle<!!id> 
  // (and recursively for all its time-like radiation products).
  // accepting only emissions which conforms to the <!id>showerVariables<!!id>
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

  static ClassDescription<BackwardEvolver> initBackwardEvolver;
  // Describe a concrete class with persistent data.

  BackwardEvolver & operator=(const BackwardEvolver &);
  //  Private and non-existent assignment operator.

  Ptr<SplittingGenerator>::pointer _splittingGenerator;
  Ptr<ForwardEvolver>::pointer _forwardEvolver;

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of BackwardEvolver.
template <>
struct BaseClassTrait<Herwig::BackwardEvolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::BackwardEvolver>
  : public ClassTraitsBase<Herwig::BackwardEvolver> {
  static string className() { return "/Herwig++/BackwardEvolver"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "BackwardEvolver.icc"

#endif /* HERWIG_BackwardEvolver_H */
