// -*- C++ -*-
#ifndef HERWIG_ForwardEvolver_H
#define HERWIG_ForwardEvolver_H
//
// This is the declaration of the ForwardEvolver class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "SplittingGenerator.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is responsible for the forward evolution of a time-like particles
 *  (and recursively to all their radiation products).
 *  It also treats the special case of the showering of a time-like 
 *  decaying particle, in which the emissions have reversed angular ordering.
 *
 * @see SplittingGenerator
 */
class ForwardEvolver: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline ForwardEvolver();
  inline ForwardEvolver(const ForwardEvolver &);
  virtual ~ForwardEvolver();

  /**
   * It does the forward evolution of the time-like input particle
   * (and recursively for all its radiation products).
   * accepting only emissions which conforms to the showerVariables
   * and soft matrix element correction pointed by meCorrectionPtr.
   * In the case that specialDecay is true then the forward evolution
   * is done with reverse angular ordering, as it should be for radiation
   * emitted by a decaying particle. 
   * If at least one emission has occurred then the method returns true
   * and all the new created ShowerParticle objects (but not the input 
   * particle) are added to the collection collecShoPar (which can
   * contain, at the beginning of the method, either the full collection
   * of ShowerParticle already created so far by the showering, 
   * or being empty: the choice is up to the caller).  
   */
  bool timeLikeShower( tEHPtr ch, 
		       const tShowerVarsPtr showerVariables, 
		       //const tMECorrectionPtr meCorrectionPtr,
		       tShowerParticlePtr particle, 
		       ShowerParticleVector & collecShoPar,
		       const bool specialDecay = false ) throw (Veto, Stop, Exception); 

private:

  bool MEVeto(tcPPtr, const Energy &, const double &);

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
  static ClassDescription<ForwardEvolver> initForwardEvolver;

  /**
   * Private and non-existent assignment operator.
   */
  ForwardEvolver & operator=(const ForwardEvolver &);

  Ptr<SplittingGenerator>::pointer _splittingGenerator;

};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ForwardEvolver.
 */
template <>
struct BaseClassTrait<Herwig::ForwardEvolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ForwardEvolver>: 
    public ClassTraitsBase<Herwig::ForwardEvolver> {
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ForwardEvolver"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwShower.so"; }
};

}

#include "ForwardEvolver.icc"

#endif /* HERWIG_ForwardEvolver_H */
