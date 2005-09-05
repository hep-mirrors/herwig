// -*- C++ -*-
#ifndef HERWIG_ScalarMesonCurrent_H
#define HERWIG_ScalarMesonCurrent_H
// This is the declaration of the ScalarMesonCurrent class.

#include "WeakDecayCurrent.h"
#include "ScalarMesonCurrent.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The weak current for the production of one (pseudo)-scalar meson.
 *
 *  In this case the current is given by
 *  \f[J^\mu = f_Pp_P^\mu,\f]
 *  where
 * - \f$f_P\f$ is the decay constant for the meson,
 * - \f$p_P\f$ is the momentum of the meson.
 *
 *  The outgoing mesons and their decay constants can be specified using the
 *  interfaces.
 *
 * @see WeakDecayCurrent.
 * 
 */
class ScalarMesonCurrent: public WeakDecayCurrent {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  ScalarMesonCurrent();

  /**
   * Copy constructor
   */
  inline ScalarMesonCurrent(const ScalarMesonCurrent &);

  /**
   * Destructor
   */
  virtual ~ScalarMesonCurrent();
  //@}

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
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector>  current(bool vertex, const int imode,
						     const int ichan, Energy & scale, 
						     const ParticleVector & decay) const;

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
   * Output the information for the database
   */
  virtual void dataBaseOutput(ofstream &) const;

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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ScalarMesonCurrent> initScalarMesonCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarMesonCurrent & operator=(const ScalarMesonCurrent &);

private:

  /**
   * the pdg code for the meson
   */
  vector<int> _id;

  /**
   * the decay constant
   */
  vector<Energy> _decay_constant;

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;

  /**
   * The inital size of the arrays
   */
  unsigned int _initsize;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ScalarMesonCurrent.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarMesonCurrent,1> {
  /** Typedef of the base class of ScalarMesonCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ScalarMesonCurrent>
  : public ClassTraitsBase<Herwig::ScalarMesonCurrent> {
  /** Return the class name. */
  static string className() { return "Herwig++::ScalarMesonCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwWeakCurrent.so"; }

};

}

#include "ScalarMesonCurrent.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMesonCurrent.tcc"
#endif

#endif /* HERWIG_ScalarMesonCurrent_H */
