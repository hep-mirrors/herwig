// -*- C++ -*-
//
// TauDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TauDecayer_H
#define HERWIG_TauDecayer_H
// This is the declaration of the TauDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"
#include "ThePEG/Helicity/LorentzSpinor.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;


/** \ingroup Decay
 *
 *  The TauDecayer class performs the decay of the \f$\tau\f$. The matrix element
 *  for \f$\tau\f$ decay can be split into a leptonic current describing the
 *  weak decay of the decay to a neutrino and a highly virtual \f$W\f$ combined
 *  with a hadronic current for the virtual \f$W\f$ decay.
 *
 *  The matrix element has the form
 *  \f[\mathcal{M} = \frac{G_F}{\sqrt{2}}\bar{u}(p_{\nu_\tau})
 *   \gamma^\mu\left(1-\gamma^5\right)u(p_{\tau}) J_\mu,\f]
 *  where
 * - \f$G_F\f$ is the Fermi constant,
 * - \f$p_{\nu_\tau}\f$ is the momentum of the \f$\tau\f$ neutrino,
 * - \f$p_{\tau}\f$ is the momentum of the \f$\tau\f$,
 * - \f$ J_\mu\f$ is the hadronic current.
 *
 *  The leptonic part of this matrix element is implemented in this class 
 *  together with a WeakCurrent member which calculates the hadronic
 *  current  \f$ J_\mu\f$. This allows a range of \f$\tau\f$ decays to be
 *  constructed via the repository using the interfaces.
 *
 * @see DecayIntegrator.
 * @see WeakCurrent
 * 
 */
class TauDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  TauDecayer() : polOpt_(false), tauMpol_(0.), tauPpol_(0.) {
    generateIntermediates(true);
  }

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
   * Output the setup information for the particle database.
   */
  void dataBaseOutput(ofstream & os,bool header) const;

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
   * Standard Init function used to initialize the interfaces.
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
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  TauDecayer & operator=(const TauDecayer &) = delete;

private:

  /**
   * mapping of the modes to the currents
   */
  vector<unsigned int> modeMap_;

  /**
   * the hadronic current
   */
  WeakCurrentPtr current_;

  /**
   * location of the weights
   */
  vector<int> wgtLoc_;

  /**
   * the maximum weight
   */
  vector<double> wgtMax_;

  /**
   *  The weights for the different channels
   */
  vector<double> weights_;

  /**
   *  The spinors for the decaying particle
   */
  mutable vector<LorentzSpinor   <SqrtEnergy> > inSpin_;

  /**
   *  Barred spinors for the deaying particle
   */
  mutable vector<LorentzSpinorBar<SqrtEnergy> > inBar_ ;

  /**
   *  Rho matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Maps for the vectors
   */
  mutable vector<unsigned int> constants_;

  /**
   *  Spins of the particles
   */
  mutable vector<PDT::Spin> iSpin_; 

  /**
   *  Option to force the polarizations of the tau leptons
   */
  bool polOpt_;

  /**
   *  Polarization for \f$\tau^-\f$
   */
  double tauMpol_;

  /**
   *  Polarization of \f$\tau^+\f$
   */
  double tauPpol_;
};

}


#endif /* HERWIG_TauDecayer_H */
