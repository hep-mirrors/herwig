// -*- C++ -*-
#ifndef HERWIG_FivePionCurrent_H
#define HERWIG_FivePionCurrent_H
//
// This is the declaration of the FivePionCurrent class.
//

#include "WeakDecayCurrent.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "FivePionCurrent.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the FivePionCurrent class.
 *
 * @see \ref FivePionCurrentInterfaces "The interfaces"
 * defined for FivePionCurrent.
 */
class FivePionCurrent: public WeakDecayCurrent {

public:

  /**
   * The default constructor.
   */
  FivePionCurrent();

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
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
   * The particles produced by the current.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector>  current(bool vertex, const int imode,
						     const int ichan,Energy & scale, 
						     const ParticleVector & decay) const;

  /**
   * Accept the decay. 
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
   *  Methods to calculate the Breit-Wigner distributions for the various
   *  mesons.
   */
  //@{
  /**
   * Breit-Wigner for the \f$\rho\f$.
   * @param scale The virtual mass
   */
  inline Complex rhoBreitWigner(Energy2 scale) const;

  /**
   * Breit-Wigner for the \f$a_1\f$.
   * @param scale The virtual mass
   */
  inline Complex a1BreitWigner(Energy2 scale) const;

  /**
   * Breit-Wigner for the \f$\omega\f$.
   * @param scale The virtual mass
   */
  inline Complex omegaBreitWigner(Energy2 scale) const;

  /**
   * Breit-Wigner for the \f$\sigma\f$.
   * @param scale The virtual mass
   */
  inline Complex sigmaBreitWigner(Energy2 scale) const;
  //@}

  /**
   *  Currents for the different channels
   */
  //@{
  /**
   * The \f$\rho\omega\f$ current
   * @param iopt Option for the inclusion of \f$\rho\f$ Breit-Wigner terms in the 
   * \f$\omega\f$ decay piece
   * @param Q The total momentum for the current
   * @param q1 The first momentum
   * @param q2 The first momentum
   * @param q3 The first momentum
   * @param q4 The first momentum
   * @param q5 The first momentum
   */
  inline LorentzPolarizationVector rhoOmegaCurrent(unsigned int iopt,
						   const Lorentz5Momentum & Q,
						   const Lorentz5Momentum & q1,
						   const Lorentz5Momentum & q2,
						   const Lorentz5Momentum & q3,
						   const Lorentz5Momentum & q4,
						   const Lorentz5Momentum & q5) const;

  /**
   *  The \f$a_1\sigma\f$ current
   * @param iopt Option for the inclusion of \f$\rho\f$ Breit-Wigner terms in the 
   * \f$a_1\f$ decay piece
   * @param Q The total momentum for the current
   * @param q1 The first momentum
   * @param q2 The first momentum
   * @param q3 The first momentum
   * @param q4 The first momentum
   * @param q5 The first momentum
   */
  inline LorentzPolarizationVector a1SigmaCurrent(unsigned int iopt,
						  const Lorentz5Momentum & Q,
						  const Lorentz5Momentum & q1,
						  const Lorentz5Momentum & q2,
						  const Lorentz5Momentum & q3,
						  const Lorentz5Momentum & q4,
						  const Lorentz5Momentum & q5) const;
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FivePionCurrent> initFivePionCurrent;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FivePionCurrent & operator=(const FivePionCurrent &);

private:

  /**
   * The masses and widths of the intermediate particles 
   */
  //@{
  /**
   * The mass of the \f$\rho\f$ for the current.
   */
  Energy _rhomass;

  /**
   * The mass of the \f$a_1\f$ for the current.
   */
  Energy _a1mass;

  /**
   * The mass of the \f$\omega\f$ for the current.
   */
  Energy _omegamass;

  /**
   * The mass of the \f$\sigma\f$ for the current.
   */
  Energy _sigmamass;

  /**
   * The width for the \f$\rho\f$.
   */
  Energy _rhowidth;

  /**
   *  The \f$a_1\f$ width
   */
  Energy _a1width;

  /**
   *  The \f$\omega\f$ width.
   */
  Energy _omegawidth;

  /**
   *  The \f$\sigma\f$ width.
   */
  Energy _sigmawidth;
  //@}

  /**
   * use local values of the particle masses
   */
  bool _localparameters;

  /**
   *  Option for the treatment of \f$\rho\f$ Breit-Wigners in \f$\omega\f$ decay
   */
  bool _rhoomega;

  /**
   *  Normalisation parameters for the different currents
   */
  //@{
  /**
   *  The \f$c\f$ parameter
   */
  Energy2 _c;

  /**
   *  The \f$c_0\f$ parameter
   */
  double _c0;

  /**
   *  The \f$f_{\omega\rho\pi}\f$ parameter
   */
  InvEnergy _fomegarhopi;

  /**
   * The \f$g_{\rho\pi\pi}\f$ parameter
   */
  double _grhopipi;

  /**
   * The \f$G_{a\rho\pi}\f$ parameter
   */
  Energy _garhopi;

  /**
   *  The \f$f_{aaf}\f$ parameter
   */
  Energy _faaf;

  /**
   *  The \f$f_{f\pi\pi}\f$ parameter
   */
  Energy _ffpipi;
  //@}

  /**
   *  Values cached to avoid unnessacary calculations
   */
  //@{
  /**
   *  Prefactor for the \f$\rho\omega\f$ current
   */
  double _preomega;

  /**
   *  Prefactor for the \f$a_1\sigma\f$ current
   */
  double _presigma;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FivePionCurrent. */
template <>
struct BaseClassTrait<Herwig::FivePionCurrent,1> {
  /** Typedef of the first base class of FivePionCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FivePionCurrent class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FivePionCurrent>
  : public ClassTraitsBase<Herwig::FivePionCurrent> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FivePionCurrent"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FivePionCurrent is implemented. It may also include several, space-separated,
   * libraries if the class FivePionCurrent depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so"; }
};

/** @endcond */

}

#include "FivePionCurrent.icc"

#endif /* HERWIG_FivePionCurrent_H */
