// -*- C++ -*-
#ifndef HERWIG_PScalarLeptonNeutrinoDecayer_H
#define HERWIG_PScalarLeptonNeutrinoDecayer_H
// This is the declaration of the PScalarLeptonNeutrinoDecayer class.

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "PScalarLeptonNeutrinoDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 * The <code>PScalarLeptonNeutrinoDecayer</code> class is designed for the decay of 
 * pseudoscalar mesons to a lepton and a neutrino. Although it can be used
 * for charged pion and kaon decays it is mainly intended for the leptonic
 * decays of bottom and charm mesons.
 *
 *  The matrix element is given by 
 * \f[\mathcal{M} = f_PG_FV_{CKM}m_l\bar{u}(p_{\ell})(1-\gamma_5)v(p_\nu),\f]
 * where
 * - \f$f_P\f$ is the pseudoscalar decay constant.
 * - \f$G_F\f$ is the Fermi constant
 * - \f$V_{CKM}\f$ is the relevant CKM matrix element
 * - \f$p_\ell\f$ is the momentum of the charged lepton
 * - \f$p_\nu\f$ is the momentum of the neutrino
 *
 * @see DecayIntegrator
 * 
 */
class PScalarLeptonNeutrinoDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline PScalarLeptonNeutrinoDecayer();

  /**
   * Copy-constructor.
   */
  inline PScalarLeptonNeutrinoDecayer(const PScalarLeptonNeutrinoDecayer &);

  /**
   * Destructor.
   */
  virtual ~PScalarLeptonNeutrinoDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks the particles against the 
   * list of allowed incoming  and outgoing mesons.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version uses PDG codes to
   * work which mode is being simulated and the generate member of the 
   * DecayIntegrator class for the phase-space  generation.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(bool vertex, const int ichan,const Particle & part,
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
  static ClassDescription<PScalarLeptonNeutrinoDecayer> initPScalarLeptonNeutrinoDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  PScalarLeptonNeutrinoDecayer & operator=(const PScalarLeptonNeutrinoDecayer &);

private:

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * the meson decay constant for a particular particle multiplied by the CKM matrix
   * element, \e i.e. \f$f_pV_{CKM}\f$
   */
  vector<Energy> _decayconstant;

  /**
   * which outgoing leptons are allowed for a particular decay
   */
  vector<int> _leptons;

  /**
   * the maximum weight for the integration of a given decay to \f$e\nu_e\f$.
   */
  vector<double> _maxweighte;

  /**
   * the maximum weight for the integration of a given decay to \f$\mu\nu_\mu\f$.
   */
  vector<double> _maxweightmu;

  /**
   * the maximum weight for the integration of a given decay to \f$\tau\nu_\tau\f$.
   */
  vector<double> _maxweighttau;

  /**
   * the fermi constant
   */
  InvEnergy2 _GF;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of PScalarLeptonNeutrinoDecayer.
 */
struct BaseClassTrait<Herwig::PScalarLeptonNeutrinoDecayer,1> {
    /** Typedef of the base class of PScalarLeptonNeutrinoDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PScalarLeptonNeutrinoDecayer>
  : public ClassTraitsBase<Herwig::PScalarLeptonNeutrinoDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::PScalarLeptonNeutrinoDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSMDecay.so"; }

};

}

#include "PScalarLeptonNeutrinoDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarLeptonNeutrinoDecayer.tcc"
#endif

#endif /* HERWIG_PScalarLeptonNeutrinoDecayer_H */
