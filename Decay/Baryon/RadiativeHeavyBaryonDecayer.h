// -*- C++ -*-
#ifndef HERWIG_RadiativeHeavyBaryonDecayer_H
#define HERWIG_RadiativeHeavyBaryonDecayer_H
//
// This is the declaration of the RadiativeHeavyBaryonDecayer class.
//

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace Herwig;

/** \ingroup Decay
 *
 * The RadiativeHeavyBaryonDecayer class is designed for the radiative decay
 * of a baryon containing a heavy quark to another baryon containing a heavy quark.
 * There are four types of transition supported
 *
 * - \f$\frac12^-\to\frac12^+\f$ \f$E1\f$ transition
 * \f[\mathcal{M} = A_{E1}\bar{B}\gamma^\mu\gamma_5B^*F^{\mu\nu}p_{0\mu}\f]
 *
 * - \f$\frac12^+\to\frac12^+\f$ \f$M1\f$ transition
 * \f[\mathcal{M} = A_{M1}\bar{B}\sigma^{\mu\nu}B^*F^{\mu\nu}\f]
 *
 * - \f$\frac32^-\to\frac12^+\f$ \f$E1\f$ transition
 * \f[\mathcal{M} = A_{E1}\bar{B}B^*_\nu F^{\mu\nu}p_{0\mu}\f]
 *
 * - \f$\frac32^+\to\frac12^+\f$ \f$M1\f$ transition
 * \f[\mathcal{M} = A_{M1}\bar{B}\gamma_\mu\gamma_5B^*_\nu F^{\mu\nu}\f]
 *
 * where \f$p_0\f$ is the momentum of the decaying baryon \f$B^*\f$ is the field for
 * the decaying baryon, \f$B\f$ is the field for the baryon produced in the decay and
 * \f$F^{\mu\nu}\f$ is the electromagnetic field strength tensor.
 *
 */
class RadiativeHeavyBaryonDecayer: public Baryon1MesonDecayerBase {

public:

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Coupling Members.
   */
  //@{
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
  virtual void threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A1,Complex& A2,Complex& A3,
					   Complex& B1,Complex& B2,Complex& B3) const;
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RadiativeHeavyBaryonDecayer & operator=(const RadiativeHeavyBaryonDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   *  The \f$M1\f$ coupling
   */
  vector<InvEnergy>  M1Coupling_;

  /**
   *  The \f$E1\f$ coupling
   */
  vector<InvEnergy2> E1Coupling_;

  /**
   * PDG code for the incoming baryons
   */
  vector<int> incoming_;

  /**
   * PDG code for the outgoing baryons
   */
  vector<int> outgoing_;

  /**
   * The type of matrix element
   */
  vector<int> modeType_;

  /**
   * max weight
   */
  vector<double> maxWeight_;
};

}

#endif /* HERWIG_RadiativeHeavyBaryonDecayer_H */
