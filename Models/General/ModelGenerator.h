// -*- C++ -*-
#ifndef HERWIG_ModelGenerator_H
#define HERWIG_ModelGenerator_H
//
// This is the declaration of the ModelGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "DecayConstructor.h"
#include "HardProcessConstructor.h"
#include "ResonantProcessConstructor.h"
#include "ModelGenerator.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to store the particles in some model and 
 * then call the appropriate function to setup the model
 *
 * @see Interfaced 
 */
class ModelGenerator: public Interfaced {

public:

  /**
   * The default constructor.
   */
  inline ModelGenerator();

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /**
   * Overloaded function from Interfaced
   */
  virtual bool preInitialize() const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}
  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ModelGenerator> initModelGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ModelGenerator & operator=(const ModelGenerator &);

private:

  /**
   * Write out the spectrum of masses and decay modes
   */
  void writeDecayModes(ofstream & ofs, tcPDPtr parent) const;

private:
  
  /**
   * Pointer to the HardProcessConstructor
   */
  HPConstructorPtr _theHPConstructor;
  
  /**
   * Pointer to DecayConstructor
   */
  DecayConstructorPtr _theDecayConstructor;

  /**
   * Vector of ParticleData pointer
   */
  PDVector _theParticles;

  /**
   * Pointer to the ResonantProcessConstructor
   */
  RPConstructorPtr _theRPConstructor;

  /**
   * The particles to create MassGenerator and WidthGenerators  
   */
  PDVector _theOffshell;
  
  /**
   * Which particles to treat as off-shell. 1 treats all particles in
   * _theParticles vector as off-shell, 0 allows selection via
   * _theOffshell vector.
   */
  int _theOffsel;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ModelGenerator. */
template <>
struct BaseClassTrait<Herwig::ModelGenerator,1> {
  /** Typedef of the first base class of ModelGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ModelGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ModelGenerator>
  : public ClassTraitsBase<Herwig::ModelGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ModelGenerator"; }
};

/** @endcond */

}

#include "ModelGenerator.icc"

#endif /* HERWIG_ModelGenerator_H */
