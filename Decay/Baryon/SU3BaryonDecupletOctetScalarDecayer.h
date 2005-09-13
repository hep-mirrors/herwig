// -*- C++ -*-
#ifndef HERWIG_SU3BaryonDecupletOctetScalarDecayer_H
#define HERWIG_SU3BaryonDecupletOctetScalarDecayer_H
//
// This is the declaration of the SU3BaryonDecupletOctetScalarDecayer class.
//
#include "Baryon1MesonDecayerBase.h"
#include "SU3BaryonDecupletOctetScalarDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>SU3BaryonDecupletOctetScalarDecayer</code> class is designed for the
 *  decay of a \f$SU(3)\f$ deculpet baryon to an octet baryon and a pseudoscalar meson
 *  from the lightest multiplet.
 *
 *  The coupling is taken to have the form
 *  \f[\frac{C}{f_\pi}\left[
 *  \bar{\Delta}^{\mu,abc}p_{\phi,\mu}\phi^i_aB^j_b\epsilon_{cij}
 * -\bar{B}^b_ip_{\phi,\mu}\phi^a_j\Delta^\mu_{abc}\epsilon^{cij}
 *  \right],\f]
 *  where \f$\Delta^\mu_{abc}\f$ is the decuplet field, \f$B^j_b\f$
 *  is the octet field and
 *  \f$\phi^i_a\f$ is the pseudoscalar field, \f$f_\pi\f$ is the pion decay constant
 *  and \f$C\f$ is the coupling for the decay.
 *
 *  This is one of a number of decayers based on \f$SU(3)\f$ symmetry which are
 *  intended for the decay of excited baryons.
 *
 * @see Baryon1MesonDecayerBase
 * @see SU3BaryonOctetOctetPhotonDecayer
 * @see SU3BaryonDecupletOctetPhotonDecayer
 * @see SU3BaryonDecupletOctetScalarDecayer
 * @see SU3BaryonSingletOctetPhotonDecayer
 * @see SU3BaryonSingletOctetScalarDecayer
 * @see SU3BaryonOctetDecupletScalarDecayer
 * @see SU3BaryonOctetOctetScalarDecayer
 * 
 */
class SU3BaryonDecupletOctetScalarDecayer: public Baryon1MesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SU3BaryonDecupletOctetScalarDecayer();

  /**
   * Copy-constructor.
   */
  inline
  SU3BaryonDecupletOctetScalarDecayer(const SU3BaryonDecupletOctetScalarDecayer &);

  /**
   * Destructor.
   */
  virtual ~SU3BaryonDecupletOctetScalarDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version tests the PDG codes against those which
   * are allowed.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   *
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version works out which
   * of the modes is required and uses the generate member of the DecayIntegrator
   * class to generate the decay.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

  /**
   *  Coupling Members.
   */
  //@{
  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A,Complex& B) const;
  //@}

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
  static ClassDescription<SU3BaryonDecupletOctetScalarDecayer>
  initSU3BaryonDecupletOctetScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SU3BaryonDecupletOctetScalarDecayer & 
  operator=(const SU3BaryonDecupletOctetScalarDecayer &);

private:

  /**
   * Set-up the modes
   */
  void setupModes(unsigned int) const;

private:

  /**
   * the coupling
   */
  double _C;

  /**
   * the relative parities of the two baryon multiplets
   */
  bool _parity;

  /**
   * the pion decay constant
   */
  Energy _fpi;

  /**
   * PDG codes for the various octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$p\f$-like member of the outgoing octet.
   */
  int _proton;

  /**
   *  The PDG code for the \f$n\f$-like member of the outgoing octet.
   */
  int _neutron;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _sigma0;

  /**
   *  The PDG code for the  \f$\Sigma^+\f$-like member of the outgoing octet.
   */
  int _sigmap;

  /**
   *  The PDG code for the  \f$\Sigma^-\f$-like member of the outgoing octet.
   */
  int _sigmam;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _lambda;

  /**
   *  The PDG code for the \f$\Xi^0\f$-like member of the outgoing octet.
   */
  int _xi0;

  /**
   *  The PDG code for the \f$\Xi^-\f$-like member of the outgoing octet.
   */
  int _xim;
  //@}

  /**
   * PDG codes for the various decuplet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$\Delta^{++}\f$-like member of the incoming decuplet.
   */
  int _deltapp;

  /**
   *  The PDG code for the \f$\Delta^{+}\f$-like member of the incoming decuplet.
   */
  int _deltap;

  /**
   *  The PDG code for the \f$\Delta^{0}\f$-like member of the incoming decuplet.
   */
  int _delta0;

  /**
   *  The PDG code for the \f$\Delta^{-}\f$-like member of the incoming decuplet.
   */
  int _deltam;

  /**
   *  The PDG code for the \f$\Sigma^{*+}\f$-like member of the incoming decuplet.
   */
  int _sigmasp;

  /**
   *  The PDG code for the \f$\Sigma^{*0}\f$-like member of the incoming decuplet.
   */
  int _sigmas0;

  /**
   *  The PDG code for the \f$\Sigma^{*-}\f$-like member of the incoming decuplet.
   */
  int _sigmasm;

  /**
   *  The PDG code for the \f$\Omega^-\f$-like member of the incoming decuplet.
   */
  int _omega;

  /**
   *  The PDG code for the \f$\Xi^{*-}\f$-like member of the incoming decuplet.
   */
  int _xism;

  /**
   *  The PDG code for the \f$\Xi^{*0}\f$-like member of the incoming decuplet.
   */
  int _xis0;
  //@}

  /**
   * PDG code for the incoming baryons
   */
  mutable vector<int> _incomingB;

  /**
   * PDG code for the outgoing baryons
   */
  mutable vector<int> _outgoingB;

  /**
   * PDG code for the outgoing mesons
   */
  mutable vector<int> _outgoingM;

  /**
   * the maximum weight for the various modes
   */
  vector<double> _maxweight;

  /**
   * The couplings for the different modes
   */
  mutable vector<InvEnergy> _prefactor;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SU3BaryonDecupletOctetScalarDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SU3BaryonDecupletOctetScalarDecayer,1> {
    /** Typedef of the base class of SU3BaryonDecupletOctetScalarDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SU3BaryonDecupletOctetScalarDecayer>
  : public ClassTraitsBase<Herwig::SU3BaryonDecupletOctetScalarDecayer> {
   /** Return the class name.*/
   static string className() { return "Herwig++::SU3BaryonDecupletOctetScalarDecayer"; }
   /**
    * Return the name of the shared library to be loaded to get
    * access to this class and every other class it uses
    * (except the base class).
    */
   static string library() { return "libHwBaryonDecay.so"; }
   
};

}

#include "SU3BaryonDecupletOctetScalarDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonDecupletOctetScalarDecayer.tcc"
#endif

#endif /* HERWIG_SU3BaryonDecupletOctetScalarDecayer_H */
