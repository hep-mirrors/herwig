// -*- C++ -*-
#ifndef Herwig_FourPionCzyzCurrent_H
#define Herwig_FourPionCzyzCurrent_H
//
// This is the declaration of the FourPionCzyzCurrent class.
//

#include "WeakDecayCurrent.h"

namespace Herwig {

using namespace ThePEG;


/** \ingroup Decay
 *
 * The FourMesonCzyzCurrent class implements the currents from Phys.Rev. D77 (2008) 114005 
 * for 4 pions
 * @see WeakDecayCurrent.
 * @see \ref FourPionCzyzCurrentInterfaces "The interfaces"
 * defined for FourPionCzyzCurrent.
 * 
 */
class FourPionCzyzCurrent: public WeakDecayCurrent {

public:

  /**
   * The default constructor.
   */
  FourPionCzyzCurrent();

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

  /**
   *   Basis current in terms of which all the others can be calculated
   */
  LorentzVector<complex<InvEnergy> > baseCurrent(Energy2 Q2,
						 const Lorentz5Momentum & Q,
						 const Lorentz5Momentum & q1,
						 const Lorentz5Momentum & q2,
						 const Lorentz5Momentum & q3,
						 const Lorentz5Momentum & q4) const;
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
  FourPionCzyzCurrent & operator=(const FourPionCzyzCurrent &);

private:
 
  /**
   *  Masses and widths of the particles
   */
  //@{
  /**
   *  Rho masses (PDG for most of current)
   */
  vector<Energy> rhoMasses_;
  /**
   *  Rho widths (PDG for most of current)
   */
  vector<Energy> rhoWidths_;
  
  /**
   *  Rho masses for the \f$F_\rho\f$ piece
   */
  vector<Energy> rhoMasses_Frho_;
  /**
   *  Rho widths for the \f$F_\rho\f$ piece
   */
  vector<Energy> rhoWidths_Frho_;
  
  /**
   *  Omega mass
   */
  Energy omegaMass_;
  /**
   *  Omega widths
   */
  Energy omegaWidth_;
  
  /**
   *  \f$f_0\f$ mass
   */
  Energy f0Mass_;
  /**
   *  \f$f_0\f$ width
   */
  Energy f0Width_;
  
  /**
   *  \f$a_1\f$ mass
   */
  Energy a1Mass_;
  /**
   *  \f$a_1\f$ width
   */
  Energy a1Width_;
  //@}

  /**
   *  Couplings in the model
   */
  //@{
  /**
   *  Coefficents for sum over \f$\rho\f$ resonances in \f$a_1\f$ term
   */
  vector<double> beta_a1_;
  
  /**
   *  Coefficents for sum over \f$\rho\f$ resonances in \f$f_0\f$ term
   */
  vector<double> beta_f0_;
  
  /**
   *  Coefficents for sum over \f$\rho\f$ resonances in \f$\omega\f$ term
   */
  vector<double> beta_omega_;
  
  /**
   *  Coefficents for sum over \f$\rho\f$ resonances in \f$B_\rho\f$ term
   */
  vector<double> beta_B_;
  
  /**
   *  Coefficents for sum over \f$\rho\f$ resonances in \f$T_\rho\f$ term
   */
  vector<double> beta_bar_;

  /**
   *   Coupling for the \f$a_1\f$ term
   */
  InvEnergy2 c_a1_;

  /**
   *   Coupling for the \f$f_0\f$ term
   */
  InvEnergy2 c_f0_;

  /**
   *   Coupling for the \f$\omega\f$ term
   */
  InvEnergy c_omega_;

  /**
   *   Coupling for the \f\rho\f$ term
   */
  double c_rho_;

  /**
   * \f$g_{\rho\pi\pi}\f$
   */
  double g_rho_pi_pi_;

  /**
   * \f$g_{\omega\pi\prho}\f$
   */
  ThePEG::Qty<std::ratio<0,1>, std::ratio<-5,1>, std::ratio<0,1> > g_omega_pi_rho_;

  /**
   * \f$g_{\rho\gamma}\f$
   */
  Energy2 g_rho_gamma_;
  //@}

  /**
   *  Pion masses
   */
  Energy mpip_, mpi0_;

};

}

#endif /* Herwig_FourPionCzyzCurrent_H */
