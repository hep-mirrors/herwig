// -*- C++ -*-
#ifndef HERWIG_StoFFFFDecayer_H
#define HERWIG_StoFFFFDecayer_H
//
// This is the declaration of the StoFFFFDecayer class.
//

#include "GeneralFourBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the StoFFFFDecayer class.
 *
 * @see \ref StoFFFFDecayerInterfaces "The interfaces"
 * defined for StoFFFFDecayer.
 */
class StoFFFFDecayer: public GeneralFourBodyDecayer {

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  StoFFFFDecayer & operator=(const StoFFFFDecayer &) = delete;

private:

  /**
   *  Signs for NO
   */
  vector<double> sign_; 
  
  /**
   *  Potential vertices for the first step
   */
  //@{
  /**
   *   VVS Vertex
   */
  vector<AbstractVVSVertexPtr> firstVVS_;

  /**
   *   VSS Vertex
   */
  vector<AbstractVSSVertexPtr> firstVSS_;

  /**
   *   SSS Vertex
   */
  vector<AbstractSSSVertexPtr> firstSSS_;

  /**
   *   FFS Vertex
   */
  vector<AbstractFFSVertexPtr> firstFFS_;
  //@}

  /**
   *  Potential vertices for the second step
   */
  //@{
  /**
   *   FFV Vertex
   */
  vector<AbstractFFVVertexPtr> secondFFV_;

  /**
   *   FFS Vertex
   */
  vector<AbstractFFSVertexPtr> secondFFS_;  
  //@}

  /**
   *  Potential vertices for the third step
   */
  //@{
  /**
   *   FFV Vertex
   */
  vector<AbstractFFVVertexPtr> thirdFFV_;

  /**
   *   FFS Vertex
   */
  vector<AbstractFFSVertexPtr> thirdFFS_;  
  //@}

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Scalar wavefunction
   */
  mutable ScalarWaveFunction swave_;

  /**
   *  Spinors for outgoing particles
   */
  mutable pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > outwave_[4];

};

}

#endif /* HERWIG_StoFFFFDecayer_H */
