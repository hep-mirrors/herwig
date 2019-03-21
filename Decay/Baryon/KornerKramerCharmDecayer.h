// -*- C++ -*-
#ifndef HERIWG_KornerKramerCharmDecayer_H
#define HERIWG_KornerKramerCharmDecayer_H
//
// This is the declaration of the KornerKramerCharmDecayer class.
//
#include "Baryon1MesonDecayerBase.h"
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

  /**
   * Default constructor.
   */
  KornerKramerCharmDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<KornerKramerCharmDecayer> initKornerKramerCharmDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  KornerKramerCharmDecayer & operator=(const KornerKramerCharmDecayer &) = delete;

private:

  /**
   * one over the number of colours 
   */
  double oneNC_;

  /**
   * Pion decay constant, \f$f_\pi\f$.
   */
  Energy fpi_;

  /**
   * Kaon decay constant, \f$f_K\f$.
   */
  Energy fk_;

  /**
   * \f$\rho\f$ decay constant, \f$f_\rho\f$.
   */
  double frho_;

  /**
   * \f$K^*\f$ decay constans, \f$f_{K^*}\f$.
   */
  double fKstar_;

  /**
   * Axial-Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to d\f$ transition.
   */
  Energy mdcplus_;

  /**
   * Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to d\f$ transition.
   */
  Energy mdcminus_;

  /**
   * Axial-Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to s\f$ transition.
   */
  Energy mscplus_;

  /**
   * Vector mass for the form factor for the factorizing diagrams
   * for the \f$c\to s\f$ transition.
   */
  Energy mscminus_;

  /**
   * Perturbative factor, \f$c_+\f$.
   */
  double cplus_;

  /**
   * Perturbative factor, \f$c_-\f$.
   */
  double cminus_;

  /**
   * \f$H_2\f$ factor for the non-factorizing diagrams.
   */
  Energy H2_;

  /**
   * \f$H_3\f$ factor for the non-factorizing diagrams.
   */
  Energy H3_;

  /**
   * SU(4) invariants for the various modes
   */
  //@{
  /**
   *  The \f$I_1\f$ invariant
   */
  vector<double> I1_;

  /**
   *  The \f$I_2\f$ invariant
   */
  vector<double> I2_;

  /**
   *  The \f$I_3\f$ invariant
   */
  vector<double> I3_;

  /**
   *  The \f$I_4\f$ invariant
   */
  vector<double> I4_;

  /**
   *  The \f$I_5\f$ invariant
   */
  vector<double> I5_;

  /**
   *  The \f$\hat{I}_3\f$ invariant
   */
  vector<double> Ihat3_;

  /**
   *  The \f$\hat{I}_4\f$ invariant
   */
  vector<double> Ihat4_;
  //@}

  /**
   * The PDG code for the incoming baryon.
   */
  vector<int> incoming_;

  /**
   * The PDG code for the outgoing baryon.
   */
  vector<int> outgoingB_;

  /**
   * The PDG code for the outgoing meson.
   */
  vector<int> outgoingM_;

  /**
   * The maximum weight.
   */
  vector<double> maxweight_;

  /**
   * The couplings for the different modes.
   */
  //@{
  /**
   * The first A coupling
   */
  vector<double> A1_;

  /**
   * The second A coupling
   */
  vector<InvEnergy> A2_;

  /**
   * The third A coupling
   */
  vector<InvEnergy2> A3_;
  /**
   * The first B coupling
   */
  vector<double> B1_;
  /**
   * The second B coupling
   */
  vector<InvEnergy> B2_;
  /**
   * The third B coupling
   */
  vector<InvEnergy2> B3_;
  //@}

  /**
   *  Initial size of the vectors
   */
  unsigned int initsize_;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

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
  static string className() { return "Herwig::KornerKramerCharmDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERIWG_KornerKramerCharmDecayer_H */
