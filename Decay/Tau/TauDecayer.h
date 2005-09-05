// -*- C++ -*-
#ifndef HERWIG_TauDecayer_H
#define HERWIG_TauDecayer_H
// This is the declaration of the TauDecayer class.

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "TauDecayer.fh"
#include "Herwig++/Decay/WeakCurrents/WeakDecayCurrent.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;


/** \ingroup Decay
 *
 *  The TauDecayer class performs the decay of the \f$\tau\f$. The matrix element
 *  for \f$\tau\f$ decay can be split into a leptonic current describing the
 *  weak decay of the decay to a neutrino and a highly virtual \f$W\f$ combined
 *  with a hadronic current for the virtual \f$W\f$ decay.
 *
 *  The matrix element has the form
 *  \f[\mathcal{M} = \frac{G_F}{\sqrt{2}}\bar{u}(p_{\nu_\tau})
 *   \gamma^\mu\left(1-\gamma^5\right)u(p_{\tau}) J_\mu,\f]
 *  where
 * - \f$G_F\f$ is the Fermi constant,
 * - \f$p_{\nu_\tau}\f$ is the momentum of the \f$\tau\f$ neutrino,
 * - \f$p_{\tau}\f$ is the momentum of the \f$\tau\f$,
 * - \f$ J_\mu\f$ is the hadronic current.
 *
 *  The leptonic part of this matrix element is implemented in this class 
 *  together with a WeakDecayCurrent member which calculates the hadronic
 *  current  \f$ J_\mu\f$. This allows a range of \f$\tau\f$ decays to be
 *  constructed via the repository using the interfaces.
 *
 * @see DecayIntegrator.
 * @see WeakDecayCurrent
 * 
 */
class TauDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline TauDecayer();

  /**
   * Copy-constructor.
   */
  inline TauDecayer(const TauDecayer &);

  /**
   * Destructor.
   */
  virtual ~TauDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This member uses the particles in the hadronic
   * current to decide if a decay is possible.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   *
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This member uses the generate
   * method from the DecayIntegrator class and the currents to perform the 
   * decay.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This method combines the leptonic current and the hadronic current to 
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
  virtual void doinit() throw(InitException);

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
  static ClassDescription<TauDecayer> initTauDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  TauDecayer & operator=(const TauDecayer &);

private:
  
  /**
   * Fermi coupling constant, \f$G_F\f$.
   */
  InvEnergy2 _GF;

  /**
   * mapping of the modes to the currents
   */
  vector<unsigned int> _modemap;

  /**
   * the hadronic current
   */
  WeakDecayCurrentPtr _current;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights and the maximums
   */
  vector<double> _wgtmax,_weights;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of TauDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::TauDecayer,1> {
    /** Typedef of the base class of TauDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TauDecayer>
  : public ClassTraitsBase<Herwig::TauDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::TauDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwWeakCurrent.so"; }

};

}

#include "TauDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TauDecayer.tcc"
#endif

#endif /* THEPEG_TauDecayer_H */
