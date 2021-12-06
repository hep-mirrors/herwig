// -*- C++ -*-
#ifndef Herwig_VectorMeson2BaryonsDecayer_H
#define Herwig_VectorMeson2BaryonsDecayer_H
//
// This is the declaration of the VectorMeson2BaryonsDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"
#include "ThePEG/Helicity/LorentzRSSpinorBar.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>VectorMeson2BaryonsDecayer</code> class is designed for the decay
 *  of a vector meson to a baryon-antibaryon pair.
 *
 *  In this case the matrix element is taken to have the form
 *  \f[\mathcal{M} = e_g \epsilon_\mu \bar{u}(p_f) \left[\gamma^\mu \right]u(p_{\bar{f}}).\f]
 *
 *  The incoming vector mesons together with their decay products and the coupling 
 *  \f$e_g\f$, \f$G_E\f$ and the phase \f$\phi_E\f$ can be specified using the interfaces for the class.
 *  The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$V\to f\bar{f}\f$ decays are specified in the default
 *  constructor.
 *
 * @see DecayIntegrator
 * @see \ref VectorMeson2BaryonsDecayerInterfaces "The interfaces"
 * defined for VectorMeson2BaryonsDecayer.
 *
 */
class VectorMeson2BaryonsDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  VectorMeson2BaryonsDecayer() {}

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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorMeson2BaryonsDecayer & operator=(const VectorMeson2BaryonsDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * \f$G_M\f$ coupling
   */
  vector<double> gm_;
  
  /**
   * \f$G_E\f$ coupling
   */
  vector<double> ge_;
  
  /**
   * Relative phase of  \f$G_E\f$/\f$G_M\f$
   */
  vector<double> phi_;

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> incoming_;

  /**
   * the PDG codes for the outgoing fermion
   */
  vector<pair<long,long> > outgoing_;

  /**
   * maximum weight for a decay
   */
  vector<double> maxweight_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors for the decaying particle
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_;

  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_;

  /**
   *  barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_;

  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzRSSpinor   <SqrtEnergy> > wave2_;

  /**
   *  barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzRSSpinorBar<SqrtEnergy> > wave2bar_;
};

}

#endif /* Herwig_VectorMeson2BaryonsDecayer_H */
