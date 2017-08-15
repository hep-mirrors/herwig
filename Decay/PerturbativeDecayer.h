// -*- C++ -*-
#ifndef Herwig_PerturbativeDecayer_H
#define Herwig_PerturbativeDecayer_H
//
// This is the declaration of the PerturbativeDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The PerturbativeDecayer class is the base class for perturbative decays in
 * Herwig and implements the functuality for the POWHEG corrections
 *
 * @see \ref PerturbativeDecayerInterfaces "The interfaces"
 * defined for PerturbativeDecayer.
 */
class PerturbativeDecayer: public DecayIntegrator {

public:
  
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}

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

protected:

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay, MEOption meopt);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PerturbativeDecayer & operator=(const PerturbativeDecayer &);

};

}

#endif /* Herwig_PerturbativeDecayer_H */
