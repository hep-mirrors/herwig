// -*- C++ -*-
#ifndef Herwig_TwoKaonCzyzCurrent_H
#define Herwig_TwoKaonCzyzCurrent_H
//
// This is the declaration of the TwoKaonCzyzCurrent class.
//

#include "WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The TwoKaonCzyzCurrent class implements the current of PRD 81 094014 for the
 * production of two kaons.
 *
 * @see \ref TwoKaonCzyzCurrentInterfaces "The interfaces"
 * defined for TwoKaonCzyzCurrent.
 */
class TwoKaonCzyzCurrent: public WeakCurrent {

public:

  /**
   * The default constructor.
   */
  TwoKaonCzyzCurrent();


  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakCurrent.
   * @param icharge   The total charge of the outgoing particles in the current.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
   * @param imode     The mode in the current being asked for.
   * @param mode      The phase space mode for the integration
   * @param iloc      The location of the of the first particle from the current in
   *                  the list of outgoing particles.
   * @param ires      The location of the first intermediate for the current.
   * @param phase     The prototype phase space channel for the integration.
   * @param upp       The maximum possible mass the particles in the current are
   *                  allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge, tcPDPtr resonance,
			  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			  unsigned int imode,PhaseSpaceModePtr mode,
			  unsigned int iloc,int ires,
			  PhaseSpaceChannel phase, Energy upp );

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
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(tcPDPtr resonance,
	  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
	  const int imode, const int ichan,Energy & scale,
	  const tPDVector & outgoing,
	  const vector<Lorentz5Momentum> & momenta,
	  DecayIntegrator::MEOption meopt) const;

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
  
  /**
   *  Calculation of the kaon form factor
   */
  Complex Fkaon(Energy2 q2,const int imode, const int ichan,
		IsoSpin::IsoSpin Itotal, tcPDPtr resonance,
		Energy ma, Energy mb) const;
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
  TwoKaonCzyzCurrent & operator=(const TwoKaonCzyzCurrent &);

private:

  /**
   * Weights for the different \f$\rho\f$ resonances in the current, \f$\alpha_k\f$.
   */
  //@{
  /**
   *  The Complex weight used in the calculation
   */
  vector<Complex> rhoWgt_;

  /**
   *  The magnitude for input
   */
  vector<double> rhoMag_;

  /**
   *  The phase for input
   */
  vector<double> rhoPhase_;
  //@}

  /**
   * Weights for the different \f$\omega\f$ resonances in the current, \f$\alpha_k\f$.
   */
  //@{
  /**
   *  The Complex weight used in the calculation
   */
  vector<Complex> omegaWgt_;

  /**
   *  The magnitude for input
   */
  vector<double> omegaMag_;

  /**
   *  The phase for input
   */
  vector<double> omegaPhase_;
  //@}

  /**
   * Weights for the different \f$\phi\f$ resonances in the current, \f$\alpha_k\f$.
   */
  //@{
  /**
   *  The Complex weight used in the calculation
   */
  vector<Complex> phiWgt_;

  /**
   *  The magnitude for input
   */
  vector<double> phiMag_;

  /**
   *  The phase for input
   */
  vector<double> phiPhase_;
  //@}
  
  /**
   * The masses of the \f$\rho\f$ resonances.
   */
  vector<Energy> rhoMasses_;

  /**
   * The widths of the \f$\rho\f$ resonances.
   */
  vector<Energy> rhoWidths_;

  /**
   * The masses of the \f$\omega\f$ resonances.
   */
  vector<Energy> omegaMasses_;

  /**
   * The widths of the \f$\omega\f$ resonances.
   */
  vector<Energy> omegaWidths_;
  
  /**
   * The masses of the \f$\phi\f$ resonances.
   */
  vector<Energy> phiMasses_;

  /**
   * The widths of the \f$\phi\f$ resonances.
   */
  vector<Energy> phiWidths_;

  /**
   *   Regge \f$\beta\f$ parameter for \f$\rho\f$ resonances
   */
  double betaRho_;

  /**
   *   Regge \f$\beta\f$ parameter for \f$\omega\f$ resonances
   */
  double betaOmega_;

  /**
   *   Regge \f$\beta\f$ parameter for \f$\phi\f$ resonances
   */
  double betaPhi_;
  
  /**
   *  Number of resonaces at which to trucated the series
   */
  unsigned int nMax_;

  /**
   *  The \f$\eta_\phi\f$ parameter
   */
  double etaPhi_;

  /**
   *  The \f$\gamma_\omega\f$ parameter
   */
  double gammaOmega_;

  /**
   *  The \f$\gamma_\phi\f$ parameter
   */
  double gammaPhi_;
  
  /**
   *   Masses of the resonances
   */
  vector<vector<Energy> > mass_;

  /**
   *   Widths of the resonances
   */
  vector<vector<Energy> > width_;

  /**
   *   Couplings of the resonaces
   */
  vector<vector<Complex> > coup_;

  /**
   * The function \f$\frac{\\hat{H}}{dq^2}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<double> dh_;

  /**
   * The function \f$\\hat{H}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<Energy2> hres_;

  /**
   * The \f$H(0)\f$ parameter  for the GS form of the
   *  Breit-Wigner
   */
  vector<Energy2> h0_;

  /**
   *  The charged pion mass
   */
  Energy mpi_;
};

}

#endif /* Herwig_TwoKaonCzyzCurrent_H */
