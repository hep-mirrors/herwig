// -*- C++ -*-
#ifndef THEPEG_VtoFFVDecayer_H
#define THEPEG_VtoFFVDecayer_H
//
// This is the declaration of the VtoFFVDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"

namespace Herwig {
  using namespace ThePEG;

/**
 * Here is the documentation of the VtoFFVDecayer class.
 *
 * @see \ref VtoFFVDecayerInterfaces "The interfaces"
 * defined for VtoFFVDecayer.
 */
class VtoFFVDecayer: public GeneralThreeBodyDecayer {

public:

  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;
  
  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of 
   * calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

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

protected:

  /**
   *   Set up the diagrams etc
   */
  virtual void setupDiagrams(bool checkKinematics);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VtoFFVDecayer & operator=(const VtoFFVDecayer &) = delete;

private:
  
  /**
   * Store the vertices for scalar intrermediate
   */
  vector<pair<AbstractVVSVertexPtr, AbstractFFSVertexPtr> > sca_;

  /**
   * Store the vertices for fermion intrermediate
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > fer_;

  /**
   * Store the vertices for vector intrermediate
   */
  vector<pair<AbstractVVVVertexPtr, AbstractFFVVertexPtr> > vec_;

  /**
   * Store the vertices for vector intrermediate
   */
  vector<pair<AbstractVVTVertexPtr, AbstractFFTVertexPtr> > ten_;

  /**
   *  Spinr density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors for the decaying particle
   */
  mutable vector<VectorWaveFunction> inVector_;

  /**
   *  Scalar wavefunction for the decay products
   */
  mutable ScalarWaveFunction swave_;

  /**
   *  Polarization vectors for the decay products
   */
  mutable vector<VectorWaveFunction> outVector_;

  /**
   *  Spinors for the decay products
   */
  mutable pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > outspin_[3];
};

}

#endif /* THEPEG_VtoFFVDecayer_H */
