// -*- C++ -*-
#ifndef Herwig_ThreePionCzyzCurrent_H
#define Herwig_ThreePionCzyzCurrent_H
//
// This is the declaration of the ThreePionCzyzCurrent class.
//

#include "WeakDecayCurrent.h"

namespace Herwig {

using namespace ThePEG;


/** \ingroup Decay
 *
 * The ThreeMesonCzyzCurrent class implements the currents from Eur.Phys.J. C47 (2006) 617-624 for 
 * \f$\pi^+\pi^-\pi^0\f$
 * @see WeakDecayCurrent.
 * @see \ref ThreePionCzyzCurrentInterfaces "The interfaces"
 * defined for ThreePionCzyzCurrent.
 * 
 */
class ThreePionCzyzCurrent: public WeakDecayCurrent {

public:

  /**
   * The default constructor.
   */
  ThreePionCzyzCurrent();

  /** @name Methods for the construction of the phase space integrator. */
  //@{ 
  /**
   * Complete the construction of the decay mode for integration.
   * This version just adds the intermediate resonances, two outgoing mesons
   * and photon.
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
   * The particles produced by the current. This just returns the two pseudoscalar
   * mesons and the photon.
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
   * @param ichan The phase-space channel the current is needed for
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode, const int ichan,Energy & scale,  
	  const ParticleVector & decay, DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the particles are the allowed mode.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
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
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ThreePionCzyzCurrent & operator=(const ThreePionCzyzCurrent &);

private:
 
  /**
   *  Masses and widths of the particles, used in the \f$I=0\f$ piece
   */
  //@{
  /**
   *  Rho masses
   */
  vector<Energy> rhoMasses_;
  /**
   *  Rho widths
   */
  vector<Energy> rhoWidths_;
  
  /**
   *  Omega masses
   */
  vector<Energy> omegaMasses_;
  /**
   *  Omega widths
   */
  vector<Energy> omegaWidths_;
  
  /**
   *  Phi mass
   */
  Energy phiMass_;
  /**
   *  Phi width
   */
  Energy phiWidth_;
  //@}

  /**
   *  Couplings in the model \f$I=0\f$, labelled A..F in paper
   */
  //@{
  vector<InvEnergy3> coup_I0_;
  //@}

  /**
   *  Masses and widths for the \f$I=1\f$ component
   */
  //@{
  /**
   *  Rho masses
   */
  vector<Energy> rhoMasses_I1_;
  /**
   *  Rho widths
   */
  vector<Energy> rhoWidths_I1_;
  
  /**
   *  Omega masses
   */
  Energy omegaMass_I1_;
  /**
   *  Omega widths
   */
  Energy omegaWidth_I1_;
  //@}
  
  /**
   *  Couplings for the the \f$I=1\f$ component
   */
  //@{
  /**
   *   The sigma parameter
   */
  double sigma_;

  /**
   *  The numerical part of \f$G_\omega\f$
   */
  InvEnergy GW_pre_;
  
  /**
   *  The full \f$G_\omega\f$
   */
  Energy GW_;

  /**
   * \f$g_{\omega\pi\pi}\f$
   */
  double g_omega_pi_pi_;
  //@}

  /**
   *  Pion mass
   */
  Energy mpip_, mpi0_;

};

}

#endif /* Herwig_ThreePionCzyzCurrent_H */
