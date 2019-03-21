// -*- C++ -*-
#ifndef HERWIG_StrongHeavyBaryonDecayer_H
#define HERWIG_StrongHeavyBaryonDecayer_H
// This is the declaration of the StrongHeavyBaryonDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>StrongHeavyBaryonDecayer</code> class implements the strong
 *  decays of charm baryons using the results of hep-ph/9904421.
 *
 * @see Baryon1MesonDecayerBase.
 * 
 */
class StrongHeavyBaryonDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  StrongHeavyBaryonDecayer();

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

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac32+0\f$ or \f$\frac32\to\frac12+0\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfThreeHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
					   Complex& A,Complex& B) const;

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

  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac32\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac32\to\frac32+0\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
						Complex& A1,Complex& A2,
						Complex& B1,Complex& B2) const;
  //@}

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
  static ClassDescription<StrongHeavyBaryonDecayer> initStrongHeavyBaryonDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  StrongHeavyBaryonDecayer & operator=(const StrongHeavyBaryonDecayer &) = delete;

private:

  /**
   * Strong coupling for the \f$\Sigma_c\to\Lambda_c\pi\f$.
   */
  InvEnergy _gsigma_clambda_cpi;

  /**
   * strong coupling for \f$\Xi^*_c\to\Xi_c\pi\f$.
   */
  InvEnergy _gxistar_cxi_cpi;

  /**
   * Strong coupling for \f$\Lambda_{c1}\to\Sigma_c\pi\f$.
   */
  double _flambda_c1sigma_cpi;

  /**
   * Strong coupling for \f$\Xi_{c1}\to\Xi'_c\pi\f$.
   */
  double _fxi_c1xi_cpi;

  /**
   * Strong coupling for \f$\Lambda_{c1}^*\to\Sigma_c\pi\f$.
   */
  InvEnergy2 _flambda_c1starsigma_cpi;

  /**
   * Strong couplng for \f$\Xi_{c1}^*\to\Xi'_c\pi\f$.
   */
  InvEnergy2 _fxi_c1starxi_cpi;

  /**
   * Strong coupling for the \f$\Sigma_b\to\Lambda_b\pi\f$.
   */
  InvEnergy _gsigma_blambda_bpi;

  /**
   * Strong coupling for \f$\Xi^*_b \to \Xi_b \pi\f$.
   */
  InvEnergy _gxistar_bxi_bpi;

  /**
   * Strong coupling for \f$\Lambda_{b1} \to \Sigma_b \pi\f$.
   */
  double _flambda_b1sigma_bpi;

  /**
   * Strong coupling for \f$\Xi_{b1}\to\Xi'_b\pi\f$.
   */
  double _fxi_b1xi_bpi;

  /**
   * Strong coupling for \f$\Lambda_{b1}^* \to \Sigma_b \pi\f$.
   */
  InvEnergy2 _flambda_b1starsigma_bpi;

  /**
   * Strong couplng for \f$\Xi_{b1}^*\to\Xi'_b\pi\f$.
   */
  InvEnergy2 _fxi_b1starxi_bpi;


  /**
   * PDG code for the incoming baryons
   */
  vector<int> _incoming;

  /**
   * PDG code for the outgoing baryons
   */
  vector<int> _outgoingB;

  /**
   * PDG code for the outgoing mesons.
   */
  vector<int> _outgoingM;

  /**
   * max weight
   */
  vector<double> _maxweight;

  /**
   * The couplings for the different modes.
   */
  vector<double> _prefactor;

  /**
   * The type of matrix element
   */
  vector<int> _modetype;

  /**
   *  The initial size of the arrays
   */
  unsigned int _initsize;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of StrongHeavyBaryonDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::StrongHeavyBaryonDecayer,1> {
    /** Typedef of the base class of StrongHeavyBaryonDecayer. */
  typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::StrongHeavyBaryonDecayer>
  : public ClassTraitsBase<Herwig::StrongHeavyBaryonDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::StrongHeavyBaryonDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_StrongHeavyBaryonDecayer_H */
