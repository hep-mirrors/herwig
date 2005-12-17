// -*- C++ -*-
#ifndef HERWIG_ShowerHandler_H
#define HERWIG_ShowerHandler_H
//
// This is the declaration of the ShowerHandler class.

#include "ThePEG/Handlers/CascadeHandler.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerParticle.h"
//#include "MECorrections.h"
#include "ShowerVariables.h"
#include "Evolver.h"

namespace Herwig {


using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is the main driver of the shower: it is responsible for 
 *  the proper handling of all other specific collaborating classes
 *  and for the storing of the produced particles in the event record.
 * 
 *  @see CascadeHandler
 *  @see InsideRangeShowerEvolver
 *  @see MECorrections
 *  @see ShowerVariables
 *  @see ShowerParticle
 */
class ShowerHandler: public ThePEG::CascadeHandler {

public:

  /**
   * Standard ctors and dtor.
   */
  inline ShowerHandler();
  inline ShowerHandler(const ShowerHandler &);
  virtual ~ShowerHandler();

public:

  /**
   * The main method which manages the all showering.
   */
  virtual void cascade();

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
  static ClassDescription<ShowerHandler> initShowerHandler;

  /**
   * Private and non-existent assignment operator.
   */
  ShowerHandler & operator=(const ShowerHandler &);

  /**
   * From the ThePEG particles entering the hard subprocess, create
   * the corresponding starting ShowerParticle objects and 
   * put them in the vector hardProcessParticles. 
   * Notice that the transformation from ThePEG ColourLine 
   * objects into ShowerColourLine ones must be properly handled.
   */
  void convertToShowerParticles(const tEHPtr ch, 
				ShowerParticleVector & hardProcessParticles);

  /**
   * It fills the positions information for all the ShowerParticle 
   * objects in _particles, at the end of the showering.
   */
  void fillPositions();

  /**
   * Print debugging information.
   */
  void debuggingInfo();

  /**
   * At the end of the Showering, transform ShowerParticle objects
   * into ThePEG particles and fill the event record with them.
   * Notice that the parent/child relationships and the 
   * transformation from ShowerColourLine objects into ThePEG
   * ColourLine ones must be properly handled.
   */
  void fillEventRecord( const tEHPtr ch );

  /**
   * Two functions to add the shower to the event record.
   * Both are recursive
   */
  void addFinalStateShower(PPtr, StepPtr);
  void addInitialStateShower(PPtr, StepPtr, bool doit=true);

  /**
   * Print the particles in the step.
   */
  void printStep(tStepPtr ptrStep, const string & title); 

  /**
   * Calculate event shape variables from a given set of particles.
   */
  void eventShape(const tShowerParticleVector & p, 
		  vector<double> & lam, vector<Vector3> & n);

  /**
   * Calculate hard ME correction.
   */
  void hardMEC(const tEHPtr ch);

  Ptr<GlobalParameters>::pointer _globalParameters; 
  //Ptr<MECorrections>::pointer _MECorrections;
  Ptr<ShowerVariables>::pointer _showerVariables;
  Ptr<Evolver>::pointer _evolver;

  /**
   *  Local storage of particles produced in shower
   */
  ShowerParticleVector _particles;   
 
};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.
 */
template <>
struct BaseClassTrait<Herwig::ShowerHandler,1> {
  typedef ThePEG::CascadeHandler NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ShowerHandler>: public ClassTraitsBase<Herwig::ShowerHandler> {
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ShowerHandler"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "ShowerHandler.icc"

#endif /* HERWIG_ShowerHandler_H */
