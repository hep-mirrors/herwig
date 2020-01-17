// -*- C++ -*-
//
// ModelGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ModelGenerator_H
#define HERWIG_ModelGenerator_H
//
// This is the declaration of the ModelGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "DecayConstructor.h"
#include "HardProcessConstructor.h"
#include "ModelGenerator.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to store the particles in some model and 
 * then call the appropriate function to setup the model
 *
 * @see \ref ModelGeneratorInterfaces "The interfaces"
 * defined for ModelGenerator.
 * @see Interfaced 
 */
class ModelGenerator: public Interfaced {

public:

  /**
   * The default constructor.
   */
  ModelGenerator() : particles_(0), offshell_(0),
		     Offsel_(0), BRnorm_(true),
		     Npoints_(10), Iorder_(1),
		     BWshape_(0), brMin_(1e-6), twoBodyOnly_(false), 
		     decayOutput_(1), minWidth_(1e-6),
		     howOffShell_(5.) {}

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
  virtual bool preInitialize() const {
    return true;
  }

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}
  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ModelGenerator & operator=(const ModelGenerator &) = delete;

private:

  /**
   * Check the decay modes a given particle type. This checks whether
   * the decay has quarks in the final state and that they can be put on
   * mass-shell during the shower.
   * @param parent The parent particle
   */
  void checkDecays(PDPtr parent);

  /**
   * Write out the spectrum of masses and decay modes
   */
  void writeDecayModes(ostream & ofs, tcPDPtr parent) const;

  /**
   * Create mass and width generators to simulate off-shell effects
   * @param p A pointer to the ParticleData object to create
   * the width and mass generators for.
   */
  void createWidthGenerator(tPDPtr p);
			   
private:
  
  /**
   * Pointer to the TwoToTwoProcessConstructor
   */
  vector<HPConstructorPtr> hardProcessConstructors_;
  
  /**
   * Pointer to DecayConstructor
   */
  DecayConstructorPtr _theDecayConstructor;

  /**
   * Vector of ParticleData pointer
   */
  PDVector particles_;

  /** @name Width and Mass Generator variables. */
  //@{
  /**
   * The particles to create MassGenerator and WidthGenerators  
   */
  PDVector offshell_;
  
  /**
   * Which particles to treat as off-shell. 1 treats all particles in
   * particles_ vector as off-shell, 0 allows selection via
   * offshell_ vector.
   */
  int Offsel_;
  
  /**
   * Whether to normalise the partial widths to BR*Total width for 
   * an on-shell particle
   */
  bool BRnorm_;

  /**
   * The number of points to include in the interpolation table
   */
  int Npoints_;
  
  /**
   * The order for the interpolation
   */
  unsigned int Iorder_;

  /**
   * The shape of the Breit-Wigner used in the mass generation
   */
  int BWshape_;

  /**
   * The minimum branching ratio to use 
   */
  double brMin_;

  /**
   *  Whether to use only two-body or all modes for running width
   */
  bool twoBodyOnly_;
  //@}

  /**
   *   Option for the outputs of the decays to a file
   */
  unsigned int decayOutput_;

  /**
   *  Minimum fraction of particle's mass width can be for off-shell
   *  treatment
   */
  double minWidth_;

  /**
   *  How much a particle is allowed to be offshell
   */
  double howOffShell_;
};

}

#endif /* HERWIG_ModelGenerator_H */
