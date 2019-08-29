// -*- C++ -*-
#ifndef Herwig_KStarKCurrent_H
#define Herwig_KStarKCurrent_H
//
// This is the declaration of the KStarKCurrent class.
//

#include "WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the KStarKCurrent class.
 *
 * @see \ref KStarKCurrentInterfaces "The interfaces"
 * defined for KStarKCurrent.
 */
class KStarKCurrent: public WeakCurrent {

public:
  
  /**
   * The default constructor.
   */
  KStarKCurrent();

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakCurrent.
   * @param icharge   The total charge of the outgoing particles in the current.
   * @param resonance If specified only include terms with this particle
   * @param flavour Information on the required flavours of the quarks
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
			  FlavourInfo flavour,
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
   * @param flavour Information on the required flavours of the quarks
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
	  FlavourInfo flavour,
	  const int imode, const int ichan,Energy & scale,
	  const tPDVector & outgoing,
	  const vector<Lorentz5Momentum> & momenta,
	  DecayIntegrator::MEOption meopt) const;

  /**
   *   Construct the SpinInfo for the decay products
   */
  virtual void constructSpinInfo(ParticleVector decay) const;

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
  
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KStarKCurrent & operator=(const KStarKCurrent &) = delete;

private:

  /**
   *  The masses of the intermediate resonances for the isoscalar piece
   */
  vector<Energy> isoScalarMasses_;
  
  /**
   *  The widths of the intermediate resonances for the isoscalar piece
   */
  vector<Energy> isoScalarWidths_;

  /**
   *  The masses of the intermediate resonances for the isovector piece
   */
  vector<Energy> isoVectorMasses_;
  
  /**
   *  The widths of the intermediate resonances for the isovector piece
   */
  vector<Energy> isoVectorWidths_;

  /**
   *  Amplitudes of the isoscalar couplings
   */
  vector<InvEnergy> isoScalarKStarAmp_;

  /**
   *  Amplitudes of the isovector couplings
   */
  vector<InvEnergy> isoVectorKStarAmp_;

  /**
   *  Phase of the isoscalar couplings
   */
  vector<double> isoScalarKStarPhase_;

  /**
   * Phase of the isovector couplings
   */
  vector<double> isoVectorKStarPhase_;

  /**
   *  Isoscalar couplings
   */
  vector<complex<InvEnergy> > isoScalarKStarCoup_;

  /**
   *  Isovector couplings
   */
  vector<complex<InvEnergy> > isoVectorKStarCoup_;

  /**
   *  The 4\f$\pi\f$ branching ratios of the resonances
   */
  vector<double> br4pi_;

  /**
   *  Branching ratio to KK*
   */
  double brKK_;

  /**
   *  Branching ratio to phieta
   */
  double brPhi_;

  /**
   *  \f$\phi\f$ mass
   */
  Energy mPhi_;

  /**
   *  \f$\eta\f$ mass
   */
  Energy mEta_;

  /**
   *  The pion mass
   */
  Energy mpi_;
};

}

#endif /* Herwig_KStarKCurrent_H */
