// -*- C++ -*-
#ifndef Herwig_WeakDalitzDecay_H
#define Herwig_WeakDalitzDecay_H
//
// This is the declaration of the WeakDalitzDecay class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Struct to contain the properties of the intermediate 
 */
struct DalitzResonance {

  /**
   *  Default constructor
   */
  DalitzResonance() {}

  /**
   *  Constructor specifiying the parameters
   */
  DalitzResonance(tPDPtr res, Energy m, Energy w,
		  unsigned int d1, unsigned int d2, unsigned int s,
		  double amp, double phi)
    : resonance(res), mass(m),width(w),
      daughter1(d1),daughter2(d2),spectator(s),
      amplitude(amp),phase(phi)
  {}

  /**
   *  Resonant particle
   */
  tPDPtr resonance;

  /**
   *  Mass of the resonance
   */
  Energy mass;

  /**
   *  Width of the resonance
   */
  Energy width;

  /**
   *  The children
   */
  unsigned int daughter1;
  unsigned int daughter2;

  /**
   *   The spectactor
   */
  unsigned int spectator;

  /**
   *  The amplitude
   */
  double amplitude;

  /**
   *  The phase
   */
  double phase;
};

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The resonance
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const DalitzResonance  & x) {
  os << x.resonance << ounit(x.mass,GeV) << ounit(x.width,GeV)
     << x.daughter1 << x.daughter2 << x.spectator
     << x.amplitude << x.phase;
  return os;
}

/** 
 * Input operator to allow the structure to be persistently read
 * @param is The input stream
 * @param x The resonance
 */
  inline PersistentIStream & operator>>(PersistentIStream & is, 
					DalitzResonance  & x) {
  is >> x.resonance >> iunit(x.mass,GeV) >> iunit(x.width,GeV)
     >> x.daughter1 >> x.daughter2 >> x.spectator
     >> x.amplitude >> x.phase;
  return is;
}
  
  
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
  WeakDalitzDecay(InvEnergy rP=5./GeV, InvEnergy rR=1.5/GeV);

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
  Complex resAmp(unsigned int i,bool gauss=false) const;

  /**
   *   Number of resonances
   */
  unsigned int nRes()  const {return resonances_.size();}
  
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
  WeakDalitzDecay & operator=(const WeakDalitzDecay &);

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
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Maximum weight for the decay
   */
  double maxWgt_;

  /**
   *  Weights for the phase-space channels
   */
  vector<double> weights_;

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
};

}

#endif /* Herwig_WeakDalitzDecay_H */
