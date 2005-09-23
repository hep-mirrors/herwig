// -*- C++ -*-
#ifndef HERWIG_NonLeptonicHyperonDecayer_H
#define HERWIG_NonLeptonicHyperonDecayer_H
//
// This is the declaration of the NonLeptonicHyperonDecayer class.
//
#include "Baryon1MesonDecayerBase.h"
#include "NonLeptonicHyperonDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This is a general class for the non-leptonic decay of hyperons. The
 *  decays are given in terms of the invariant amplitudes
 *
 *  \f[\bar{u}_{B_j} \left\{A+B\gamma_5\right\}u_{B_i}\f]
 *
 *  where \f$B_j\f$ is the outgoing baryon and \f$B_i\f$ is the incoming baryon.
 *
 *  The default amplitudes are taken from the fit in hep-ph/9902351, 
 *  N.B. due to the sign conventions in hep-ph/9902351 the B amplitudes
 *  have the opposite sign.
 *
 * @see Baryon1MesonDecayerBase
 * 
 */
class NonLeptonicHyperonDecayer: public Baryon1MesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline NonLeptonicHyperonDecayer();

  /**
   * Copy-constructor.
   */
  inline NonLeptonicHyperonDecayer(const NonLeptonicHyperonDecayer &);

  /**
   * Destructor.
   */
  virtual ~NonLeptonicHyperonDecayer();
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
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
				      Complex& A,Complex& B) const;

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

protected:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<NonLeptonicHyperonDecayer> initNonLeptonicHyperonDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  NonLeptonicHyperonDecayer & operator=(const NonLeptonicHyperonDecayer &);

private:

  /**
   * PDG code for the incoming baryon.
   */
  vector<int> _incomingB;

  /**
   * PDG code for the outgoing baryon.
   */
  vector<int> _outgoingB;

  /**
   * PDG code for the outgoing meson
   */
  vector<int> _outgoingM;

  /**
   * The \f$A\f$ coefficient.
   */
  vector<double> _A;

  /**
   * The \f$B\f$ coefficient.
   */
  vector<double> _B;

  /**
   * the maximum weights for the decays
   */
  vector<double> _maxweight;

  /**
   *  initial size fo the vectors
   */
  unsigned int _initsize;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of NonLeptonicHyperonDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::NonLeptonicHyperonDecayer,1> {
    /** Typedef of the base class of NonLeptonicHyperonDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::NonLeptonicHyperonDecayer>
  : public ClassTraitsBase<Herwig::NonLeptonicHyperonDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig++::NonLeptonicHyperonDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

}

#include "NonLeptonicHyperonDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "NonLeptonicHyperonDecayer.tcc"
#endif

#endif /* HERWIG_NonLeptonicHyperonDecayer_H */
