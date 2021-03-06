// -*- C++ -*-
#ifndef HERWIG_BaryonFactorizedDecayer_H
#define HERWIG_BaryonFactorizedDecayer_H
//
// This is the declaration of the BaryonFactorizedDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"
#include "Herwig/Decay/FormFactors/BaryonFormFactor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/Helicity/LorentzRSSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The BaryonFactorizedDecayer class is designed to combine the form factor
 * for a weak baryon transition and a weak decay current to produce a decayer.
 * It is mainly based on the results of PRD56, 2799.
 *
 * @see BaryonFactorizedDecayer
 * @see BaryonFormFactor
 * @see WeakCurrent
 */
class BaryonFactorizedDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  BaryonFactorizedDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

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

protected:

  /**
   * Matrix element for \f$\frac12\to\frac12\f$.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfHalf(const int ichan,
		  const Particle & part,
		  const tPDVector & outgoing,
		  const vector<Lorentz5Momentum> & momenta,
		  MEOption meopt) const;

  /**
   * Matrix element for \f$\frac12\to\frac32\f$.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfThreeHalf(const int ichan,
		       const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const;
  
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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Find duplicate modes in the list of particles
   * @param imode The mode we are studying
   * @param incoming The incoming particles for the different modes
   * @param outgoing The outgoing particles for the different modes
   * @param loc The location of the duplicate mode
   * @param cc  If the duplicate is the charge conjugate
   */
  void findModes(unsigned int imode,
		 tPDVector & incoming,
		 vector<tPDVector> & outgoing,
		 vector<unsigned int> & loc,vector<bool> & cc);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BaryonFactorizedDecayer & operator=(const BaryonFactorizedDecayer &) = delete;

private:

  /**
   *  The weak decay current
   */
  WeakCurrentPtr _current;

  /**
   * The baryon form factor
   */
  BaryonFormFactorPtr _form;

  /**
   *  The perturbative coefficients
   */
  //@{
  /**
   *  The perturbative \f$a_1\f$ coefficient for b decays.
   */
  double _a1b;

  /**
   *  The perturbative \f$a_2\f$ coefficient for b decays.
   */
  double _a2b;

  /**
   *  The perturbative \f$a_1\f$ coefficient for c decays.
   */
  double _a1c;

  /**
   *  The perturbative \f$a_2\f$ coefficient for c decays.
   */
  double _a2c;
  //@}

  /**
   * Mapping of the modes to the currents
   */
  vector<vector<unsigned int> > _currentmap;

  /**
   * Mapping of the modes to the form factors
   */
  vector<vector<unsigned int> > _formmap;

  /**
   * The CKM factors
   */
  vector<vector <Complex> > _factCKM;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights
   */
  vector<double> _wgtmax;

  /**
   *  weights for the different channels
   */
  vector<double> _weights;

  /**
   * Pointer to the CKM object.
   */
  Ptr<StandardCKM>::pointer _theCKM;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *   Spin-\f$\frac12\f$ spinors
   */
  mutable vector<LorentzSpinor<SqrtEnergy> > _inHalf;

  /**
   *   Spin-\f$\frac12\f$ barred spinors
   */
  mutable vector<LorentzSpinorBar<SqrtEnergy> > _inHalfBar;

  /**
   *   Spin-\f$\frac32\f$ spinors
   */
  mutable vector<LorentzRSSpinor<SqrtEnergy> > _inThreeHalf;

  /**
   *   Spin-\f$\frac32\f$ barred spinors
   */
  mutable vector<LorentzRSSpinorBar<SqrtEnergy> > _inThreeHalfBar;
};

}

#endif /* HERWIG_BaryonFactorizedDecayer_H */
