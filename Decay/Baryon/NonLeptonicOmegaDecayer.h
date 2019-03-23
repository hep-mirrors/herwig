// -*- C++ -*-
#ifndef HERWIG_NonLeptonicOmegaDecayer_H
#define HERWIG_NonLeptonicOmegaDecayer_H
// This is the declaration of the NonLeptonicOmegaDecayer class.

#include "Baryon1MesonDecayerBase.h"

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

  /**
   * Default constructor.
   */
  NonLeptonicOmegaDecayer();

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
   * Private and non-existent assignment operator.
   */
  NonLeptonicOmegaDecayer & operator=(const NonLeptonicOmegaDecayer &) = delete;

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
  double _cbstar;

  /**
   * The \f$s_c\f$ coupling of the \f$R\f$ multiplet to the decuplet.
   */
  double _sc;

  /**
   * The \f$\mathcal{C}\f$ coupling of the decuplet to the ground state baryons.
   */
  double _c;

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
  Energy _mlambda;

  /**
   * The mass of the \f$\Xi\f$.
   */
  Energy _mxi;

  /**
   * The mass of the \f$\Omega\f$.
   */
  Energy _momega;

  /**
   * The mass of the \f$\Xi^*\f$.
   */
  Energy _mxistar;

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
  Energy _mkp;

  /**
   * The mass of the \f$K^0\f$.
   */
  Energy _mk0;

  /**
   * The mass of the \f$B^*\f$ resonance, this is the \f$\frac12^+\f$ multiplet.
   */
  Energy _mbstar;

  /**
   * The mass of the \f$R\f$ resonance, this is the \f$\frac12^-\f$ multiplet.
   */
  Energy _mr;

  /**
   * use local values for the masses for the couplings
   */
  bool _localmasses;

  /**
   * The PDG code for the incoming baryon.
   */
  long _incomingB;

  /**
   * The PDG code for the outgoing baryon.
   */
  vector<long> _outgoingB;

  /**
   *  The PDG code for the outgoing meson.
   */
  vector<long> _outgoingM;

  /**
   * The \f$A\f$ coefficient for the decays.
   */
  vector<InvEnergy> _a;

  /**
   * The \f$B\f$ coefficient for the decays.
   */
  vector<InvEnergy> _b;

  /**
   * The maximum weights for the decays.
   */
  vector<double> _maxweight;
};

}


#endif /* HERWIG_NonLeptonicOmegaDecayer_H */
