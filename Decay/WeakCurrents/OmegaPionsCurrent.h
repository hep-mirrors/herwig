// -*- C++ -*-
#ifndef HERWIG_OmegaPionsCurrent_H
#define HERWIG_OmegaPionsCurrent_H
//
// This is the declaration of the OmegaPionsCurrent class.
//

#include "WeakDecayCurrent.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "OmegaPionsCurrent.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The OmegaPionsCurrent class is designed to implement the results of
 * PRD55, 1436 1997, hep-ph/9607354 and hep-ph/0004097 for the 
 * weak current of \f$\omega\f$ with one, two or three pions.
 *
 * @see WeakDecayCurrent
 * @see \ref OmegaPionsCurrentInterfaces "The interfaces"
 * defined for OmegaPionsCurrent.
 */
class OmegaPionsCurrent: public WeakDecayCurrent {

public:

  /**
   * The default constructor.
   */
  inline OmegaPionsCurrent();

public:

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
   * Hadronic current. 
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector> current(bool vertex, const int imode,
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
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<OmegaPionsCurrent> initOmegaPionsCurrent;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OmegaPionsCurrent & operator=(const OmegaPionsCurrent &);

private:

  /**
   *  Breit-Wigner for the \f$a_1\f$.
   * @param iopt Include the numerator (iopt=0) or not (iopt=1)
   * @param q2 The scale
   */
  inline Complex a1BreitWigner(unsigned int iopt, Energy q2) const;

  /**
   *  Breit-Wigner for the \f$\rho\f$.
   * @param iopt Include the numerator (iopt=0) or not (iopt=1)
   * @param q2 The scale
   */
  inline Complex rhoBreitWigner(unsigned int iopt, Energy q2) const;

  /**
   *  The running width for the \f$\rho\f$
   * @param q The scale.
   */
  inline Energy rhoWidth(Energy q) const;

  /**
   * The running width for the \f$\rho\f$
   * @param q The scale.
   */
  inline Energy a1Width(Energy q) const;

  /**
   * The \f$A(q^2,p^2)\f$ function for the \f$a_1\rho\pi\f$ vertex 
   * @param q2 The mass of the \f$a_1\f$.
   * @param p2 The mass of the \f$\rho\f$.
   */
  inline Energy Acoupling(Energy2 q2, Energy2 p2) const; 

  /**
   *  The \f$f^(12)_{\alpha\lambda}\f$ function contracted with a vector for
   *  the alphas index and multiplied by the \f$\rho\f$ BreitWigner factor
   * @param alpha The vector to contract with the first index
   * @param q2 The scale for the current
   * @param q The total momentum for the current
   * @param p1 The momentum of the first particle
   * @param p2 The momentum of the second particle
   */
  LorentzPolarizationVector ffunction(LorentzPolarizationVector & alpha,
				      Energy2 q2,
				      const Lorentz5Momentum & q,
				      const Lorentz5Momentum & p1,
				      const Lorentz5Momentum & p2) const;
private:

  /**
   *  The pion decay constant \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   *  The coupling \f$g\f$.
   */
  double _gfact;

  /**
   *  The \f$\rho\f$ mass
   */
  Energy _mrho;

  /**
   *  The \f$a_1\f$ mass
   */
  Energy _ma1;

  /**
   *  Use local values of the \f$\rho\f$ and \f$a_1\f$ mass,
   */
  bool _localparameters;

  /**
   *   The \f$a_1\f$ coupling
   */
  double _fa;

  /**
   *   The \f$c\f$ parameter
   */
  double _cfact;

  /**
   *  the pion mass
   */
  Energy _mpi;

  /**
   *  the \f$B\f$ coupling for the \f$a_1\rho\pi\f$ vertex 
   */
  InvEnergy _Bcoup;

  /**
   *  the \f$D\f$ coupling for the \f$a_1\rho\pi\f$ vertex 
   */
  InvEnergy _Dcoup;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of OmegaPionsCurrent. */
template <>
struct BaseClassTrait<Herwig::OmegaPionsCurrent,1> {
  /** Typedef of the first base class of OmegaPionsCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the OmegaPionsCurrent class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::OmegaPionsCurrent>
  : public ClassTraitsBase<Herwig::OmegaPionsCurrent> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::OmegaPionsCurrent"; }
  /**
   * The name of a file containing the dynamic library where the class
   * OmegaPionsCurrent is implemented. It may also include several, space-separated,
   * libraries if the class OmegaPionsCurrent depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so"; }
};

/** @endcond */

}

#include "OmegaPionsCurrent.icc"

#endif /* HERWIG_OmegaPionsCurrent_H */
