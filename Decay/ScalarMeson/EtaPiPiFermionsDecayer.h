// -*- C++ -*-
//
// EtaPiPiFermionsDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EtaPiPiFermionsDecayer_H
#define HERWIG_EtaPiPiFermionsDecayer_H
// This is the declaration of the EtaPiPiFermionsDecayer class.

#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "Herwig/Decay/FormFactors/OmnesFunction.fh"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>EtaPiPiFermionsDecayer</code> class implements the decay of
 * the \f$\eta\f$ or \f$\eta'\f$ to \f$\pi^+\pi^-\gamma\f$ using either 
 * a VMD type model or a model using either the theoretical or experimental 
 * form of the Omnes function taken from hep-ph/0112150.
 *
 * The matrix element is given by
 * \f[\mathcal{M} = B(s_{+-},s_{+\gamma},s_{-\gamma})\epsilon^{\mu\nu\alpha\beta}
 *      \epsilon^*_{\mu}p_{+\nu}p_{-\alpha}p_{\gamma\beta}\f]
 * where \f$p_{+,-}\f$ are the momenta of the positively and negatively charged pions,
 * \f$p_{\gamma}\f$ is the momentum of the photon and \f$s_{ij} = (p_i+p_j)^2\f$.
 *
 *  The different models take
 *
 *  \f[B(s_{+-},s_{+\gamma},s_{-\gamma}) = 
 *  B_0\left(1+\frac32\frac{s_{+-}}{M^2_\rho-s_{+-}-iM_\rho\Gamma_\rho(s_{+-})}\right)\f]
 *  where \f$M_\rho\f$ and \f$\Gamma_\rho\f$ are the mass and running width 
 *  of the \f$\rho\f$
 *  respectively for the VMD model.
 *
 *  For the Omnes function case we take
 *
 *  \f[B(s_{+-},s_{+\gamma},s_{-\gamma}) = 
 *  B_0\left(1-c+c\frac{1+as_{+-}}{D_1(s_{+-})}\right)\f]
 *  either the experimental or analytic form of the Omnes function \f$D_1(s_{+-})\f$
 *  taken from hep-ph/0112150 can be used.
 *
 *  The coefficient \f$B_0\f$ is given in hep-ph/0112150. We use the values from this
 *  paper and use their default choice \f$c=1\f$, \f$a=\frac1{2M_\rho}\f$.
 *
 * @see DecayIntegrator
 * 
 */
class EtaPiPiFermionsDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  EtaPiPiFermionsDecayer();
  
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

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * Private and non-existent assignment operator.
   */
  EtaPiPiFermionsDecayer & operator=(const EtaPiPiFermionsDecayer &) = delete;

private:

  /**
   * the pion decay constant, \f$F_\pi\f$.
   */
  Energy fPi_;

  /**
   * the PDG code for the incoming particle
   */
  vector<int> incoming_;

  /**
   *  The PDG code of the leptons
   */
  vector<int> lepton_;

  /**
   * Coupling for the decay, \f$B_0\f$.
   */
  vector<double> coupling_;

  /**
   * The maximum weight
   */
  vector<double> maxWeight_;

  /**
   * The option for the energy dependence of the prefactor
   */
  vector<int> option_;

  /**
   * The constants for the omnes function form.
   */
  InvEnergy2 aConst_;

  /**
   * The constants for the Omnes function form.
   */
  double cConst_;

  /**
   * The \f$\rho\f$ mass
   */
  Energy mRho_;

  /**
   * The \f$\rho\f$ width
   */
  Energy rhoWidth_;

  /**
   * Constant for the running \f$rho\f$ width.
   */
  double rhoConst_;

  /**
   * The \f$m_\pi\f$.
   */
  Energy mPi_;

  /**
   * Use local values of the parameters.
   */
  bool localParameters_;

  /**
   *  Object calculating the Omnes function
   */
  OmnesFunctionPtr omnesFunction_;

  /**
   *  Spin densit matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Spinors for the fermions
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_;

  /**
   *  Barred spinors for the fermions 
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_;
};
}

#endif /* HERWIG_EtaPiPiFermionsDecayer_H */
