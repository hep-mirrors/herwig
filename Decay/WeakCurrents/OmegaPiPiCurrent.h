// -*- C++ -*-
#ifndef Herwig_OmegaPiPiCurrent_H
#define Herwig_OmegaPiPiCurrent_H
//
// This is the declaration of the OmegaPiPiCurrent class.
//

#include "WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the OmegaPiPiCurrent class.
 *
 * @see \ref OmegaPiPiCurrentInterfaces "The interfaces"
 * defined for OmegaPiPiCurrent.
 */
class OmegaPiPiCurrent: public WeakCurrent {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  OmegaPiPiCurrent() {
  // /**
  //  *   Parameters for the \f$\omega(1650)\f$
  //  */
  // //@{
  // /**
  //  *   Mass of the resonance
  //  */
  // mRes_;
  
  // /**
  //  *   Width of the resonance
  //  */
  // wRes_;

  // /**
  //  *   Coupling of the resonance
  //  */
  // gRes_;
  // //@}

  // /**
  //  *   Parameters for the \f$f_0\f$ resonances
  //  */
  // //@{
  // /**
  //  *  Mass of the \f$\sigma\f$
  //  */
  // mSigma_;
  
  // /**
  //  *  Width of the \f$\sigma\f$
  //  */
  // wSigma_;
  
  // /**
  //  *  Mass of the \f$f_0(980)\f$
  //  */
  // mf0_;

  // /**
  //  *  \f$f_0\f$ coupling to \f$\pi\pi\f$
  //  */
  // gPiPi_;

  // /**
  //  *  \f$f_0\f$ coupling to KK
  //  */
  // gKK_;

  // /**
  //  *   Sigma coupling
  //  */
  // gSigma_;
  
  // /**
  //  *   f_0 coupling
  //  */
  // gf0_;
  };
  //@}


public:

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
   *   Construct the SpinInfo for the decay products
   */
  virtual void constructSpinInfo(ParticleVector decay) const;

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
   * Initialize this object after the setup phase before saving an
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
  OmegaPiPiCurrent & operator=(const OmegaPiPiCurrent &) = delete;

private:

  /**
   *   Parameters for the \f$\omega(1650)\f$
   */
  //@{
  /**
   *   Mass of the resonance
   */
  Energy mRes_;
  
  /**
   *   Width of the resonance
   */
  Energy wRes_;

  /**
   *   Coupling of the resonance
   */
  double gRes_;
  //@}

  /**
   *   Parameters for the \f$f_0\f$ resonances
   */
  //@{
  /**
   *  Mass of the \f$\sigma\f$
   */
  Energy mSigma_;
  
  /**
   *  Width of the \f$\sigma\f$
   */
  Energy wSigma_;
  
  /**
   *  Mass of the \f$f_0(980)\f$
   */
  Energy mf0_;

  /**
   *  \f$f_0\f$ coupling to \f$\pi\pi\f$
   */
  double gPiPi_;

  /**
   *  \f$f_0\f$ coupling to KK
   */
  double gKK_;

  /**
   *   Sigma coupling
   */
  double gSigma_;
  
  /**
   *   f_0 coupling
   */
  double gf0_;
  //@}
};

}

#endif /* Herwig_OmegaPiPiCurrent_H */
