// -*- C++ -*-
#ifndef HERWIG_BaryonFactorizedDecayer_H
#define HERWIG_BaryonFactorizedDecayer_H
//
// This is the declaration of the BaryonFactorizedDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/WeakCurrents/WeakDecayCurrent.h"
#include "Herwig++/Decay/FormFactors/BaryonFormFactor.h"
#include "BaryonFactorizedDecayer.fh"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The BaryonFactorizedDecayer class is designed to combine the form factor
 * for a weak baryon transition and a weak decay current to produce a decayer.
 * It is mainly based on the results of PRD56, 2799.
 *
 * @see BaryonFactorizedDecayer
 * @see BaryonFormFactor
 * @see WeakDecayCurrent
 */
class BaryonFactorizedDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline BaryonFactorizedDecayer();

  /**
   * The copy constructor.
   */
  inline BaryonFactorizedDecayer(const BaryonFactorizedDecayer &);

  /**
   * The destructor.
   */
  virtual ~BaryonFactorizedDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This method combines the form factor and the weka current to 
   * calculate the matrix element.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		      const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database.
   */
  void dataBaseOutput(ofstream &) const;

protected:

  /**
   * Matrix element for \f$\frac12\to\frac12\f$.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param inpart The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfHalf(bool vertex, const int ichan,const Particle & inpart,
		  const ParticleVector & decay) const;

  /**
   * Matrix element for \f$\frac12\to\frac32\f$.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param inpart The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfThreeHalf(bool vertex, const int ichan,const Particle & inpart,
		       const ParticleVector & decay) const;

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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
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
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Find duplicate modes in the list of particles
   * @param imode The mode we are studying
   * @param particles The external particles for the different modes
   * @param loc The location of the duplicate mode
   * @param cc  If the duplicate is the charge conjugate
   */
  void findModes(unsigned int imode,vector<PDVector> & particles,
		 vector<unsigned int> & loc,vector<bool> & cc);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BaryonFactorizedDecayer> initBaryonFactorizedDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BaryonFactorizedDecayer & operator=(const BaryonFactorizedDecayer &);

private:

  /**
   *  The weak decay current
   */
  WeakDecayCurrentPtr _current;

  /**
   * The baryon form factor
   */
  BaryonFormFactorPtr _form;

  /**
   *  The Fermi constant, \f$G_F\f$
   */
  InvEnergy2 _GF;

  /**
   *  The perturbative coefficients
   */
  //@{
  /**
   *  The perturbative \f$a_1\f$ coefficient for b decays.
   */
  double _a1b;

  /**
   *  The perturbative \f$a_2\f$ coefficient for b decays.
   */
  double _a2b;

  /**
   *  The perturbative \f$a_1\f$ coefficient for c decays.
   */
  double _a1c;

  /**
   *  The perturbative \f$a_2\f$ coefficient for c decays.
   */
  double _a2c;
  //@}

  /**
   * Mapping of the modes to the currents
   */
  vector<vector<unsigned int> > _currentmap;

  /**
   * Mapping of the modes to the form factors
   */
  vector<vector<unsigned int> > _formmap;

  /**
   * The CKM factors
   */
  vector<vector <Complex> > _CKMfact;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights and the weights
   */
  vector<double> _wgtmax,_weights;

  /**
   * Pointer to the CKM object.
   */
  Ptr<StandardCKM>::pointer _theCKM;
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of BaryonFactorizedDecayer. */
template <>
struct BaseClassTrait<Herwig::BaryonFactorizedDecayer,1> {
  /** Typedef of the first base class of BaryonFactorizedDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BaryonFactorizedDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BaryonFactorizedDecayer>
  : public ClassTraitsBase<Herwig::BaryonFactorizedDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::BaryonFactorizedDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BaryonFactorizedDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwBaryonDecay.so"; }
};

}

#include "BaryonFactorizedDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonFactorizedDecayer.tcc"
#endif

#endif /* HERWIG_BaryonFactorizedDecayer_H */
