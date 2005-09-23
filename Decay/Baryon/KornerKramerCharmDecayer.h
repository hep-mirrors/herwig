// -*- C++ -*-
#ifndef HERIWG_KornerKramerCharmDecayer_H
#define HERIWG_KornerKramerCharmDecayer_H
//
// This is the declaration of the KornerKramerCharmDecayer class.
//
#include "Baryon1MesonDecayerBase.h"
#include "KornerKramerCharmDecayer.fh"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>KornerKramerCharmDecayer</code> class implements the model of 
 *  Z.Phys.C55,659 (1992) for the non-leptonic decay of charm baryons.
 *  The couplings of the model are calculated at initialisation and stored. These
 *  couplings are then returned when requested using the coupling members of the 
 *  base class.
 *
 * @see Baryon1MesonDecayerBase.
 * 
 */
class KornerKramerCharmDecayer: public Baryon1MesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  KornerKramerCharmDecayer();

  /**
   * Copy-constructor.
   */
  inline KornerKramerCharmDecayer(const KornerKramerCharmDecayer &);

  /**
   * Destructor.
   */
  virtual ~KornerKramerCharmDecayer();
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

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a vector.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void halfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
				      Complex& A1,Complex& A2,
				      Complex& B1,Complex& B2) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a vector.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param A3 The coupling \f$A_3\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   * @param B3 The coupling \f$B_3\f$ described above.
   */
  virtual void halfThreeHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A1,Complex& A2,Complex& A3,
					   Complex& B1,Complex& B2,Complex& B3) const;
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
  static ClassDescription<KornerKramerCharmDecayer> initKornerKramerCharmDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  KornerKramerCharmDecayer & operator=(const KornerKramerCharmDecayer &);

private:

  /**
   * The Fermi constant, \f$G_F\f$.
   */
  InvEnergy2 _GF;

  /**
   * one over the number of colours 
   */
  double _oneNC;

  /**
   * Pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * Kaon decay constant, \f$f_K\f$.
   */
  Energy _FK;

  /**
   * \f$\rho\f$ decay constant, \f$f_\rho\f$.
   */
  double _frho;

  /**
   * \f$K^*\f$ decay constans, \f$f_{K^*}\f$.
   */
  double _fKstar;

  /**
   * Axial-Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to d\f$ transition.
   */
  Energy _mdcplus;

  /**
   * Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to d\f$ transition.
   */
  Energy _mdcminus;

  /**
   * Axial-Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to s\f$ transition.
   */
  Energy _mscplus;

  /**
   * Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to s\f$ transition.
   */
  Energy _mscminus;

  /**
   * Perturbative factor, \f$c_+\f$.
   */
  double _cplus;

  /**
   * Perturbative factor, \f$c_-\f$.
   */
  double _cminus;

  /**
   * \f$H_2\f$ factor for the non-factorizing diagrams.
   */
  Energy _H2;

  /**
   * \f$H_3\f$ factor for the non-factorizing diagrams.
   */
  Energy _H3;

  /**
   * SU(4) invariants for the various modes
   */
  //@{
  /**
   *  The \f$I_1\f$ invariant
   */
  vector<double> _I1;

  /**
   *  The \f$I_2\f$ invariant
   */
  vector<double> _I2;

  /**
   *  The \f$I_3\f$ invariant
   */
  vector<double> _I3;

  /**
   *  The \f$I_4\f$ invariant
   */
  vector<double> _I4;

  /**
   *  The \f$I_5\f$ invariant
   */
  vector<double> _I5;

  /**
   *  The \f$\hat{I}_3\f$ invariant
   */
  vector<double> _Ihat3;

  /**
   *  The \f$\hat{I}_4\f$ invariant
   */
  vector<double> _Ihat4;
  //@}

  /**
   * The PDG code for the incoming baryon.
   */
  vector<int> _incoming;

  /**
   * The PDG code for the outgoing baryon.
   */
  vector<int> _outgoingB;

  /**
   * The PDG code for the outgoing meson.
   */
  vector<int> _outgoingM;

  /**
   * The maximum weight.
   */
  vector<double> _maxweight;

  /**
   * The couplings for the different modes.
   */
  //@{
  /**
   * The first A coupling
   */
  vector<double> _A1;

  /**
   * The second A coupling
   */
  vector<InvEnergy> _A2;

  /**
   * The third A coupling
   */
  vector<InvEnergy2> _A3;
  /**
   * The first B coupling
   */
  vector<double> _B1;
  /**
   * The second B coupling
   */
  vector<InvEnergy> _B2;
  /**
   * The third B coupling
   */
  vector<InvEnergy2> _B3;
  //@}

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of KornerKramerCharmDecayer.
 */
  template <>
  struct BaseClassTrait<Herwig::KornerKramerCharmDecayer,1> {
    /** Typedef of the base class of KornerKramerCharmDecayer. */
    typedef Herwig::Baryon1MesonDecayerBase NthBase;
  };

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::KornerKramerCharmDecayer>
  : public ClassTraitsBase<Herwig::KornerKramerCharmDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::KornerKramerCharmDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

}

#include "KornerKramerCharmDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KornerKramerCharmDecayer.tcc"
#endif

#endif /* HERIWG_KornerKramerCharmDecayer_H */
