// -*- C++ -*-
#ifndef Herwig_f1RhoPiPiDecayer_H
#define Herwig_f1RhoPiPiDecayer_H
//
// This is the declaration of the f1RhoPiPiDecayer class.
//

#include "DecayIntegrator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the f1RhoPiPiDecayer class.
 *
 * @see \ref f1RhoPiPiDecayerInterfaces "The interfaces"
 * defined for f1RhoPiPiDecayer.
 */
class f1RhoPiPiDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  f1RhoPiPiDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;
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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  f1RhoPiPiDecayer & operator=(const f1RhoPiPiDecayer &);

private:

  /**
   *  rho-pi-pi coupling
   */
  double gRhoPiPi_;

  /**
   *  a_1-rho-pi
   */
  Energy ga1RhoPi_;

  /**
   *  f1-a1-pi
   */
  InvEnergy gf1a1Pi_;

  /**
   *  a1 mass
   */
  Energy ma1_;

  /**
   * a1 width
   */
  Energy ga1_;
  

  /**
   *  Weights for the different decay modes
   */
  vector<double> maxWeight_;

  /**
   *  Spin Density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_[2];
  
};

}

#endif /* Herwig_f1RhoPiPiDecayer_H */
