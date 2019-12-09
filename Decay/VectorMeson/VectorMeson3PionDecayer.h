// -*- C++ -*-
//
// VectorMeson3PionDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VectorMeson3PionDecayer_H
#define HERWIG_VectorMeson3PionDecayer_H
// This is the declaration of the VectorMeson3PionDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The VectorMeson3PionDecayer class is designed to perform the
 *  decay of an \f$I=0\f$ meson to three pions via rho mesons including the
 *  option of higher rho resonaces and a constant term. It is mainly intended for
 *  the decays:
 *   - \f$\omega \to \pi^+\pi^-\pi^0\f$;
 *   - \f$\phi   \to \pi^+\pi^-\pi^0\f$.
 *
 *  The default for the \f$\omega\f$ is to only include the
 *  contributions of the \f$\rho(770)\f$
 *  without a constant term, whereas for the \f$\phi\f$ the parameters
 *  from hep-ex/0303016 (KLOE) which includes the \f$\rho(770)\f$ 
 *  and a constant term to represent the effects
 *  of the higher rho resonances is used.
 * (The KLOE paper also included a omega contribution but this is assumed to be a 
 *  non-resonant contribution to the \f$e^+e^-\f$ cross section.)
 *
 *  The form of the matrix element is taken to be
 *
 *  \f[\mathcal{M} = \frac{g}{M^2_{\rho}}\epsilon^{\mu\alpha\beta\gamma}
 *                      \epsilon_{0\mu}
 *                      p_{+\alpha}p_{-\beta}p_{0\gamma}
 *      \left[a_de^{i\phi_d} +\sum_k a_{\rho_k}e^{i\phi_{\rho_k}}
 *   \left\{
 *  \frac{M^2_{\rho_k^+}}{m^2_{0+}-M^2_{\rho_k^+}
 *           +im_{0+}\Gamma_{\rho_k^+}(m^2_{0+})}
 *  +\frac{M^2_{\rho_k^-}}{m^2_{0-}-M^2_{\rho_k^-}
 *           +im_{0-}\Gamma_{\rho_k^-}(m^2_{0-})}
 *  +\frac{M^2_{\rho_k^0}}{m^2_{+-}-M^2_{\rho_k^0}
 *           +im_{+-}\Gamma_{\rho_k^0}(m^2_{+-})}\right\}
 *  \right],\f]
 *  where \f$\epsilon_0\f$ is the polarization vector of the decaying meson,
 *  \f$p_{+,-,0}\f$ are the momenta of the positively charged, negatively charged
 *  and neutral pions respectively, \f$m^2_{ij}=(p_i+p_j)^2\f$, \f$M_{\rho_k}\f$ is
 *  the mass of the \f$k\f$th \f$\rho\f$ resonance, \f$g\f$ is the overall coupling and 
 * \f[\Gamma_{\rho_k}(m^2) = \Gamma_{\rho_k}\left(\frac{p_\pi(m^2)}{p_\pi(M_{\rho_k}^2)}\right)^2
 *    \left(\frac{M_{\rho_k}^2}{m^2}\right) \f]
 *  where \f$p_\pi\f$ is the pion momentum in the \f$\rho\f$ rest frame
 *  and \f$\Gamma_{\rho_k}\f$ is the width. In practice the couplings of the different
 *  terms are measured relative to the coupling of the lightest \f$\rho\f$ multiplet
 *   which is assumed to be one.
 *
 *
 *  To allow the easy addition of further modes the parameters for additional modes
 *  can be set. The following must be specified for each mode
 *
 *   - Incoming       - the PDG code for the incoming particle
 *   - Coupling       - the overall coupling for the decay, \f$g\f$.
 *   - DirectCoupling - the relative coupling for the direct term, \f$a_d\f$.
 *   - DirectPhase    - the phase of the coupling for the direct term, \f$\phi_d\f$.
 *   - Rho2Coupling   - the relative coupling for the second rho multiplet, 
 *     \f$a_{\rho_2}\f$.
 *   - Rho2Phase      - the phase of the coupling for the second rho multiplet, 
 *     \f$\phi_{\rho_2}\f$.
 *   - Rho3Coupling   - the relative coupling for the third rho multiplet, 
 *     \f$a_{\rho_3}\f$.
 *   - Rho3Phase      - the phase of the coupling for the third rho multiplet, 
 *     \f$\phi_{\rho_3}\f$.
 *   - MaximumWeight  - the maximum weight for the integration of the channel
 *   - Rho1Weight     - the weight for the first rho multiplet in the 
 *                      multichannel integration
 *   - Rho2Weight     - the weight for the second rho multiplet in the
 *                      multichannel integration
 *   - Rho3Weight     - the weight for the third rho multiplet in the
 *                      multichannel integration
 *   - Rho1Mass       - mass  of the first  rho multiplet, \f$M_{\rho_1}\f$.
 *   - Rho2Mass       - mass  of the second rho multiplet, \f$M_{\rho_2}\f$.
 *   - Rho3Mass       - mass  of the third  rho multiplet, \f$M_{\rho_3}\f$.
 *   - Rho1Width      - width of the first  rho multiplet, \f$\Gamma_{\rho_1}\f$.
 *   - Rho2Width      - width of the second rho multiplet, \f$\Gamma_{\rho_2}\f$.
 *   - Rho3Width      - width of the third  rho multiplet, \f$\Gamma_{\rho_3}\f$.
 *
 * @see DecayIntegrator
 * @see \ref VectorMeson3PionDecayerInterfaces "The interfaces"
 * defined for VectorMeson3PionDecayer.
 * 
 *  \author Peter Richardson
 */
class VectorMeson3PionDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  VectorMeson3PionDecayer();

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(const int imode, const Energy2 q2,
					const  Energy2 s3, const Energy2 s2, const 
					Energy2 s1, const Energy m1,
					const Energy m2, const Energy m3) const;
  
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

private:

  /**
   * Private and non-existent assignment operator.
   */
  VectorMeson3PionDecayer & operator=(const VectorMeson3PionDecayer &) = delete;

private:

  /**
   * vector storing the decaying particles for the different modes
   */
  vector<double> _incoming;

  /**
   * the overall coupling for the decay, \f$g\f$.
   */
  vector<InvEnergy> _coupling;

  /**
   * relative coupling for the direct term, \f$a_d\f$.
   */
  vector<double> _directcoupling;

  /**
   * relative phase for the direct term, \f$\phi_d\f$.
   */
  vector<double> _directphase;

  /**
   * relative coupling for the second rho multiplet, \f$a_{\rho_2}\f$.
   */
  vector<double> _rho2coupling;

  /**
   * relative phase for the second rho multiplet, \f$\phi_{\rho_2}\f$.
   */
  vector<double> _rho2phase;

  /**
   * relative coupling for the second rho multiplet, \f$a_{\rho_3}\f$.
   */
  vector<double> _rho3coupling;

  /**
   * relative phase for the second rho multiplet, \f$\phi_{\rho_3}\f$.
   */
  vector<double> _rho3phase;

  /**
   * maximum weight for the integration of the channel
   */
  vector<double> _maxwgt;

  /**
   * weight for the first  rho multiplet in the integration
   */
  vector<double> _rho1wgt;

  /**
   * weight for the second rho multiplet in the integration
   */
  vector<double> _rho2wgt;

  /**
   * weight for the third  rho multiplet in the integration
   */
  vector<double> _rho3wgt;

  /**
   * mass of the first rho multiplet, \f$M_{\rho_1}\f$.
   */
  vector<Energy> _rho1mass;

  /**
   * mass of the  second rho multiplet, \f$M_{\rho_2}\f$.
   */
  vector<Energy> _rho2mass;

  /**
   * mass of the third rho multiplet, \f$M_{\rho_3}\f$.
   */
  vector<Energy> _rho3mass;

  /**
   * width of the first rho multiplet, \f$\Gamma_{\rho_1}\f$.
   */
  vector<Energy> _rho1width;

  /**
   * width of the  second rho multiplet, \f$\Gamma_{\rho_2}\f$.
   */
  vector<Energy> _rho2width;

  /**
   * width of the third rho multiplet, \f$\Gamma_{\rho_3}\f$.
   */
  vector<Energy> _rho3width;

  /**
   * use the default parameters for the rho masses and widths
   */
  vector<bool> _defaultmass;

  /**
   * Constants for the running widths
   */
  //@{
  /**
   *  For the neutral \f$\rho\f$
   */
  vector<vector <double> > _rho0const;

  /**
   *  For the charged \f$\rho\f$
   */
  vector<vector <double> > _rhocconst;
  //@}

  /**
   * rho mass parameters
   */
  vector<vector<Energy> > _rhomass;

  /**
   * rho mass parameters
   */
  vector<vector<Energy2> > _rhomass2;

  /**
   * couplings as complex numbers
   */
  vector<vector <complex<InvEnergy2> > > _ccoupling;

  /**
   * The charge pion mass
   */
  Energy _mpic;

  /**
   *  The neutral pion mass
   */
  Energy _mpi0;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

  /**
   *  Storage of polarization tensors to try and increase
   *  speed
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors;
  
  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix _rho;
};

}


#endif /* HERWIG_VectorMeson3PionDecayer_H */
