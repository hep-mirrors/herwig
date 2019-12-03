// -*- C++ -*-
//
// TwoPionPhotonCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoPionPhotonCurrent_H
#define HERWIG_TwoPionPhotonCurrent_H
//
// This is the declaration of the TwoPionPhotonCurrent class.
//
#include "WeakDecayCurrent.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This class implements the decay current for \f$\pi^\pm\pi^0 \gamma\f$ via
 *  an intermediate \f$\omega\f$. It inherits from the <code>WeakDecayCurrent</code>
 *  class and implements the hadronic current.
 *
 *  The model is based on the one used in TAUOLA, Comput.Phys.Commun.76:361-380,1993.
 *  The current is given by
 * \f[J^\mu  = e T \left\{
 *     \epsilon^\mu\left[ m^2_\pi      p_1\cdot p_3
 *                       -p_2\cdot p_3(p_2\cdot p_1-p_1\cdot p_3)\right]
 *     -p_2^\mu\left[p_2\cdot\epsilon p_1\cdot p_3-p_1\cdot\epsilon p_2\cdot p_3\right]
 *     +p_3^\mu\left[\epsilon\cdot p_2-\epsilon\cdot p_1(m^2_\pi+p_2\cdotp_3)\right]
 *\right\}\f]
 *  where
 * - \f$p_1\f$ is the momentum of the charged pion
 * - \f$p_2\f$ is the momentum of the neutral pion
 * - \f$p_3\f$ is the momentum of the photon
 * - \f$\epsilon\f$ is the polarization of the photon
 * - \f$e\f$ is the electric charge of the positron
 * and the normaliztion factor is
 * \f[T = F(q^2)F(0)\frac1{\sqrt{2}B_\omega(s_2)}\f]
 * and
 * \f[F(s) = \sqrt{2}F_\rho g_{\rho\omega\pi}\sum_k\sigma_k B_{\rho_k}(s)\f]
 * where
 * - \f$B_\omega(s)=\frac1{m^2_\omega-s-im_\omega\Gamma_\omega}\f$ is the Breit-Wigner for the \f$\omega\f$.
 * - \f$m_\omega\f$ is the mass of the \f$\omega\f$.
 * - \f$\Gamma_\omega \f$ is the width of the \f$\omega\f$.
 * - \f$F_\rho\f$ is the coupling for the conversion of the \f$\rho\f$ to a photon.
 * - \f$g_{\rho\omega\pi}\f$ is the coupling of \f$\rho\f$, \f$\omega\f$, \f$\pi\f$.
 * - \f$m_{\rho_k}\f$ is the mass of the \f$k\f$th \f$\rho\f$ resonance
 * - \f$B_{\rho_k}(s)=\frac1{m^2_{\rho_k}-s-im_{\rho_k}\Gamma_{\rho_k}}\f$ is the 
 *    Breit-Wigner for  \f${\rho_k}\f$.
 * - \f$m_{\rho_k}\f$ is the mass of the \f${\rho_k}\f$.
 * - \f$\Gamma_{\rho_k} \f$ is the width of the \f${\rho_k}\f$.
 *
 * @see WeakDecayCurrent
 * 
 *  \author Peter Richardson
 *
 */
class TwoPionPhotonCurrent: public WeakDecayCurrent {

public:

  /**
   * Default constructor
   */
  TwoPionPhotonCurrent();

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

public:


  /** @name Methods for the construction of the phase space integrator. */
  //@{ 
  /**
   * Complete the construction of the decay mode for integration.
   * This version just adds the meson as the daughter of the last
   * resonance in the phase space channel.
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp);

  /**
   * The particles produced by the current. This just returns the pseudoscalar
   * meson.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode,const int ichan,Energy & scale, 
	  const ParticleVector & decay,DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the meson against the list
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Checks the meson against the list
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  TwoPionPhotonCurrent & operator=(const TwoPionPhotonCurrent &) = delete;

private:
  
  /**
   * Calculate the \f$F(q^2)\f$ function at a given scale
   * @param q2 The scale \f$q^2\f$.
   * @return The value of the function. 
   */
  complex<InvEnergy> FFunction(Energy2 q2) const {
    complex<InvEnergy2> output(ZERO);
    for(unsigned int ix=0, N=_resweights.size(); ix<N && ix <3;++ix) {
      output -= _resweights[ix]*BreitWigner(q2,ix);
    }
    return output*_grho*_grhoomegapi*sqrt(2.);
  }

  /**
   * Fixed width Breit wigner
   * @param q2 The scame \f$q^2\f$
   * @param ires The resonance required (0,1,3) are the \f$\rho\f$'s and 10 is the 
   * \f$\omega\f$.
   * @return The breit wigner
   */
  complex<InvEnergy2> BreitWigner(Energy2 q2,unsigned int ires) const {
    static const Complex ii(0.,1.);
    complex<Energy2> denom;
    if(ires<_rhomasses.size()) {
      denom = q2-_rhomasses[ires]*_rhomasses[ires]+ii*_rhomasses[ires]*_rhowidths[ires];
    }
    else if(ires==10) {
      denom = q2-_omegamass*_omegamass+ii*_omegamass*_omegawidth;
    }
    else assert(false);
    return 1./denom;
  }
  
private:
  
  /**
   * Coupling of the rho to the photon, \f$F_\rho\f$.
   */
  Energy2 _grho;
  
  /**
   * Coupling of the rho to the omega and a pion, \f$g_{\rho\omega\pi}\f$.
   */
  InvEnergy _grhoomegapi;
  
  /**
   * Weights of the different rho resonances in the current
   */
  vector<double> _resweights;
  
  /**
   * Use local parameters for the rho resonances rather than from the particle data
   * objects
   */
  bool _rhoparameters;

  /**
   * Masses of the \f$\rho\f$ resonances
   */
  vector<Energy> _rhomasses; 

  /**
   * Widths of the \f$\rho\f$ resonances
   */
  vector<Energy> _rhowidths; 

  /**
   * use local parameters for the omega rather than from the particle data objects
   */
  bool _omegaparameters;

  /**
   * The \f$\omega\f$ mass.
   */
  Energy _omegamass;

  /**
   * The \f$\omega\f$ width.
   */
  Energy _omegawidth;

  /**
   * Mass for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intmass;

  /**
   * Width for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intwidth; 

};

}


#endif /* HERWIG_TwoPionPhotonCurrent_H */
