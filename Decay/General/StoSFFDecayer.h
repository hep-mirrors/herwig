// -*- C++ -*-
#ifndef THEPEG_StoSFFDecayer_H
#define THEPEG_StoSFFDecayer_H
//
// This is the declaration of the StoSFFDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSSVertex.h"

namespace Herwig {
  using namespace ThePEG;

/**
 * The StoSFFDecayer class provides the general matrix element for
 * scalar decays into another scalar and a fermion-antifermion pair.
 *
 * @see \ref StoSFFDecayerInterfaces "The interfaces"
 * defined for StoSFFDecayer.
 */
class StoSFFDecayer: public GeneralThreeBodyDecayer {

public:

  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;
  
  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
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
  StoSFFDecayer & operator=(const StoSFFDecayer &) = delete;

private:
  
  /**
   * Store the vertices for scalar intermediate
   */
  vector<pair<AbstractSSSVertexPtr, AbstractFFSVertexPtr> > sca_;

  /**
   * Store the vertices for spin-\f$\frac12\f$ fermion intermediate
   */
  vector<pair<AbstractFFSVertexPtr, AbstractFFSVertexPtr> > fer_;

  /**
   * Store the vertices for spin-\f$\frac32\f$  fermion intermediate
   */
  vector<pair<AbstractRFSVertexPtr, AbstractRFSVertexPtr> > RSfer_;

  /**
   * Store the vertices for vector intermediate
   */
  vector<pair<AbstractVSSVertexPtr, AbstractFFVVertexPtr> > vec_;

  /**
   * Store the vertices for tensor intermediate
   */
  vector<pair<AbstractSSTVertexPtr, AbstractFFTVertexPtr> > ten_;

  /**
   * Store the vertices for four point diagrams
   */
  vector<AbstractFFSSVertexPtr> four_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Scalar wavefunction
   */
  mutable ScalarWaveFunction swave_;

  /**
   *  Spinor wavefunctions
   */
  mutable pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > outspin_[3];
};

}

#endif /* THEPEG_StoSFFDecayer_H */
