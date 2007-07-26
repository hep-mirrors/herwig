// -*- C++ -*-
#ifndef HERWIG_NonLeptonicOmegaDecayer_H
#define HERWIG_NonLeptonicOmegaDecayer_H
// This is the declaration of the NonLeptonicOmegaDecayer class.

#include "Baryon1MesonDecayerBase.h"
#include "NonLeptonicOmegaDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>NonLeptonicOmegaDecayer</code> class is designed for the non-leptonic
 *  weak decay of the Omega to a baryon from the lightest \f$SU(3)\f$ octet and a 
 *  pseudoscalar meson. The results are taken from hep-ph/9905398.
 *
 *  due to problems with the size of the d-wave term and recent measurements giving
 * the opposite sign for the \f$\alpha\f$ parameter we have set this term to zero.
 *
 * @see Baryon1MesonDecayerBase.
 * 
 */
class NonLeptonicOmegaDecayer: public Baryon1MesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline NonLeptonicOmegaDecayer();

  /**
   * Copy-constructor.
   */
  inline NonLeptonicOmegaDecayer(const NonLeptonicOmegaDecayer &);

  /**
   * Destructor.
   */
  virtual ~NonLeptonicOmegaDecayer();
  //@}

public:

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool & cc,const DecayMode & dm) const;

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
  static ClassDescription<NonLeptonicOmegaDecayer> initNonLeptonicOmegaDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  NonLeptonicOmegaDecayer & operator=(const NonLeptonicOmegaDecayer &);

private:

  /**
   * The \f$d^*\f$ coupling for the \f$B^*\f$ multiplet, this is \f$d^* /M_{B^*}\f$
   * from the paper.
   */
  double _dstar;

  /**
   * The \f$f^*\f$ coupling for the \f$B^*\f$ multiplet, this is \f$f^* /M_{B^*}\f$
   * from the paper.
   */
  double _fstar;

  /**
   * The \f$\omega_d\f$ coupling for the \f$R\f$ multiplet, this is \f$\omega_d/M_R\f$
   * from the paper.
   */
  double _omegad;

  /**
   * The \f$\omega_f\f$ coupling for the \f$R\f$ multiplet, this is \f$\omega_f/M_R\f$
   * from the paper.
   */
  double _omegaf;

  /**
   * The \f$\mathcal{C}_{B^*}\f$ coupling of the \f$B^*\f$ multiplet to the decuplet.
   */
  double _CBstar;

  /**
   * The \f$s_c\f$ coupling of the \f$R\f$ multiplet to the decuplet.
   */
  double _sc;

  /**
   * The \f$\mathcal{C}\f$ coupling of the decuplet to the ground state baryons.
   */
  double _C;

  /**
   * The pion decay constant \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * The \f$h_\pi\f$ pion self-coupling.
   */
  double _hpi;

  /**
   * The \f$h_c\f$ coupling.
   */
  Energy _hc;

  /**
   * The \f$d\f$ coupling for the ground-state baryon multiplet.
   */
  Energy _d;

  /**
   * The \f$f\f$ coupling for the ground-state baryon multiplet.
   */
  Energy _f;

  /**
   * The mass of the \f$\Lambda^0\f$.
   */
  Energy _Mlambda;

  /**
   * The mass of the \f$\Xi\f$.
   */
  Energy _Mxi;

  /**
   * The mass of the \f$\Omega\f$.
   */
  Energy _Momega;

  /**
   * The mass of the \f$\Xi^*\f$.
   */
  Energy _MXistar;

  /**
   * The mass of the \f$\pi^+\f$.
   */
  Energy _mpip;

  /**
   * The mass of the \f$\pi^0\f$.
   */
  Energy _mpi0;

  /**
   * The mass of the \f$K^+\f$.
   */
  Energy _MKp;

  /**
   * The mass of the \f$K^0\f$.
   */
  Energy _MK0;

  /**
   * The mass of the \f$B^*\f$ resonance, this is the \f$\frac12^+\f$ multiplet.
   */
  Energy _MBstar;

  /**
   * The mass of the \f$R\f$ resonance, this is the \f$\frac12^-\f$ multiplet.
   */
  Energy _MR;

  /**
   * use local values for the masses for the couplings
   */
  bool _localmasses;

  /**
   * The PDG code for the incoming baryon.
   */
  int _incomingB;

  /**
   * The PDG code for the outgoing baryon.
   */
  vector<int> _outgoingB;

  /**
   *  The PDG code for the outgoing meson.
   */
  vector<int> _outgoingM;

  /**
   * The \f$A\f$ coefficient for the decays.
   */
  vector<InvEnergy> _A;

  /**
   * The \f$B\f$ coefficient for the decays.
   */
  vector<InvEnergy> _B;

  /**
   * The maximum weights for the decays.
   */
  vector<double> _maxweight;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of NonLeptonicOmegaDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::NonLeptonicOmegaDecayer,1> {
    /** Typedef of the base class of NonLeptonicOmegaDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::NonLeptonicOmegaDecayer>
  : public ClassTraitsBase<Herwig::NonLeptonicOmegaDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::NonLeptonicOmegaDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#include "NonLeptonicOmegaDecayer.icc"

#endif /* HERWIG_NonLeptonicOmegaDecayer_H */
