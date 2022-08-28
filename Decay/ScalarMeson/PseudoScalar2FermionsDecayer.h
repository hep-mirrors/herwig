// -*- C++ -*-
#ifndef Herwig_PseudoScalar2FermionsDecayer_H
#define Herwig_PseudoScalar2FermionsDecayer_H
//
// This is the declaration of the PseudoScalar2FermionsDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The PseudoScalar2FermionsDecayer class is designed for the decay of a pseudoscalar meson to
 * two fermions, in reality this is always the decay of the \f$\eta_c\f$ to a baryon
 * antibaryon pair.
 *
 * @see \ref PseudoScalar2FermionsDecayerInterfaces "The interfaces"
 * defined for PseudoScalar2FermionsDecayer.
 */
class PseudoScalar2FermionsDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  PseudoScalar2FermionsDecayer();

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
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 2.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  PseudoScalar2FermionsDecayer & operator=(const PseudoScalar2FermionsDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * coupling for a decay
   */
  vector<double> coupling_;

  /**
   * the PDG codes for the incoming particles
   */
  vector<long> incoming_;

  /**
   * the PDG codes for the outgoing fermion-antifermion
   */
  vector<pair<long,long>> outgoing_;

  /**
   * maximum weight for a decay
   */
  vector<double> maxweight_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_;

  /**
   *  barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_;

};

}

#endif /* Herwig_PseudoScalar2FermionsDecayer_H */
