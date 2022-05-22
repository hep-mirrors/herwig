// -*- C++ -*-
#ifndef Herwig_WeakDalitzDecay_H
#define Herwig_WeakDalitzDecay_H
//
// This is the declaration of the WeakDalitzDecay class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "DalitzResonance.h"

namespace Herwig {
using namespace ThePEG;
  
  
/**
 * The WeakDalitzDecay class provides a base class for the implementation
 * of weak three-body decays of bottom and charm mesons
 *
 * @see \ref WeakDalitzDecayInterfaces "The interfaces"
 * defined for WeakDalitzDecay.
 */
class WeakDalitzDecay: public DecayIntegrator {
  
public:

  /**
   * The default constructor.
   */
  WeakDalitzDecay(InvEnergy rP=5./GeV, InvEnergy rR=1.5/GeV, bool useResonanceMass=false);

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
   *  Add a new resonance
   */
  void addResonance(const DalitzResonance & R) {resonances_.push_back(R);}

  /**
   *  Set up the phase-space
   */
  void createMode(tPDPtr in, tPDVector out);

  /**
   *  Calculate the amplitude
   */
  virtual Complex amplitude(int ichan) const = 0;

  /**
   * Calculate the amplitude for the ith resonance
   */
  Complex resAmp(unsigned int i) const;

  /**
   *  Access to the resonances
   */
  const vector<DalitzResonance> & resonances() const {return resonances_;}

  /**
   *  Access to the invariants
   */
  const Energy & mInv(unsigned int i,unsigned int j) const {
    return m2_[i][j];
  }

  /**
   *  Masses of the children
   */
  const Energy & mOut(unsigned int i) const {
    return mOut_[i];
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
  WeakDalitzDecay & operator=(const WeakDalitzDecay &) = delete;

private:

  /**
   *   The radii for the Blatt-Weisskopf form factors
   */
  //@{
  /**
   *   For the decaying particles
   */
  InvEnergy rParent_;

  /**
   *  For the intermediate resonances
   */
  InvEnergy rResonance_;
  //@}

  /**
   *  Vector containing the intermediate resonances
   */
  vector<DalitzResonance> resonances_;

  /**
   *  Choice of the mass to use in the denominator of the expressions
   */
  bool useResonanceMass_;
  
private:
  
  /**
   *   Parameters for the phase-space sampling
   */
  //@{
  /**
   *  Maximum weight for the decay
   */
  double maxWgt_;

  /**
   *  Weights for the phase-space channels
   */
  vector<double> weights_;
  //@}
  
private :

  /**
   *   Storage of the kinematics
   */
  //@{
  /**
   *   Mass of the parent
   */
  mutable Energy mD_;

  /**
   *   Masses of the children
   */
  mutable Energy mOut_[3];

  /**
   *   Masses of the children
   */
  mutable Energy m2_[3][3];
  //@}

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;
};

}

#endif /* Herwig_WeakDalitzDecay_H */
