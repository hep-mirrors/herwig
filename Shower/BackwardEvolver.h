// -*- C++ -*-
#ifndef HERWIG_BackwardShowerEvolver_H
#define HERWIG_BackwardShowerEvolver_H
//
// This is the declaration of the BackwardShowerEvolver class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/PartialCollisionHandler.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerConfig.h"
#include "SplittingGenerator.h"
#include "ForwardEvolver.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * 
 *  This class is responsible for the backward evolution of a 
 *  space-like particles (and recursively to all their time-like 
 *  radiation products).
 *
 *  @see SplittingGenerator
 *  @see ForwardShowerEvolver
 */
class BackwardEvolver: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline BackwardEvolver();
  inline BackwardEvolver(const BackwardEvolver &);
  virtual ~BackwardEvolver();

  /**
   * It does the backward evolution of the space-like input particle 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables
   * and soft matrix element correction pointed by meCorrectionPtr.
   * The ParticleCollisionHandler object is needed to access the PDFs.
   * If at least one emission has occurred then the method returns true
   * and all the new created ShowerParticle objects (but not the input
   * particle) are added to the collection collecShoPar (which can
   * contain, at the beginning of the method, either the full collection
   * of ShowerParticle already created so far by the showering, 
   * or being empty: the choice is up to the caller).
   */
  bool spaceLikeShower( tPartCollHdlPtr ch, 
		        const tShowerVarsPtr showerVariables, 
		        //const tMECorrectionPtr meCorrectionPtr,
		        tShowerParticlePtr particle, 
			ShowerParticleVector &allShowerParticles) 
    throw (Veto, Stop, Exception);

private:

  /**
   * This routine is used to generate the ShowerKinematics object in the 
   * forced splitting.
   */
  ShoKinPtr forcedSplitting(const ShowerParticle &, Energy, Energy);

  /**
   * This routine sets all the properties of the new particles from the 
   * splitting: it fixes the hadron parent/children relations due to the 
   * splitting and the colour information.
   */
  void createBranching(ShowerParticlePtr, ShowerParticlePtr, 
		       ShowerParticlePtr, Energy, 
		       ShowerIndex::InteractionType);

  /**
   * This routine sets the colour connections for a backwards branching.
   */
  void setColour(ShowerParticlePtr&, ShowerParticlePtr&, ShowerParticlePtr&);
		 
public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<BackwardEvolver> initBackwardEvolver;

  /**
   * Private and non-existent assignment operator.
   */
  BackwardEvolver & operator=(const BackwardEvolver &);

  Ptr<SplittingGenerator>::pointer _splittingGenerator;
  Ptr<ForwardEvolver>::pointer _forwardEvolver;

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of BackwardEvolver.
 */
template <>
struct BaseClassTrait<Herwig::BackwardEvolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::BackwardEvolver>
  : public ClassTraitsBase<Herwig::BackwardEvolver> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/BackwardEvolver"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwShower.so"; }

};

}

#include "BackwardEvolver.icc"

#endif /* HERWIG_BackwardEvolver_H */
