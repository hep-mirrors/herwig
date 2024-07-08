// -*- C++ -*-
#ifndef HERWIG_KornerKurodaFormFactor_H
#define HERWIG_KornerKurodaFormFactor_H
//
// This is the declaration of the KornerKurodaFormFactor class.
//

#include "BaryonFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The KornerKurodaFormFactor class implements the simple model of PRD 16, 2165, 1977 for the
 * form factors of the baryon octet/
 *
 * @see \ref KornerKurodaFormFactorInterfaces "The interfaces"
 * defined for KornerKurodaFormFactor.
 */
class KornerKurodaFormFactor: public BaryonFormFactor {

public:

  /**
   * The default constructor.
   */
  KornerKurodaFormFactor();

  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac12\f$ baryon.  
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming baryon.
   * @param id1 The PDG code of the outgoing baryon.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param f1v The form factor \f$F^V_1\f$.
   * @param f2v The form factor \f$F^V_2\f$.
   * @param f3v The form factor \f$F^V_3\f$.
   * @param f1a The form factor \f$F^A_1\f$.
   * @param f2a The form factor \f$F^A_2\f$.
   * @param f3a The form factor \f$F^A_3\f$.
   * @param virt Whether the q2 is space or timelike
   */
  virtual void SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
					  Energy m0, Energy m1,
					  Complex & f1v,Complex & f2v,Complex & f3v,
					  Complex & f1a,Complex & f2a,Complex & f3a,
					  FlavourInfo flavour,
					  Virtuality virt=SpaceLike);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
  
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KornerKurodaFormFactor & operator=(const KornerKurodaFormFactor &) = delete;

private:

  /**
   *  Whether or not to include the nucleon modes
   */
  bool includeNucleon_;

  /**
   * Masses for the form-factors
   */
  //@{
  /**
   *  Mass of the \f$\rho\f$ meson
   */
  Energy mRho_;

  /**
   *  Mass of the \f$\omega\f$ meson
   */
  Energy mOmega_;

  /**
   *  Mass of the \f$\phi\f$ meson
   */
  Energy mPhi_;
  //@}

  /**
   *   Regge slope
   */
  InvEnergy2 aPrime_;

  /**
   *  Couplings for the form factors
   */
  //@{
  /**
   *  Coupling \f$C_1^\rho\f$
   */
  vector<double> c1Rho_;
  
  /**
   *  Coupling \f$C_1^\omega\f$
   */
  vector<double> c1Omega_;
  
  /**
   *  Coupling \f$C_1^\phi\f$
   */
  vector<double> c1Phi_;
  
  /**
   *  Coupling \f$C_1^\rho+C_2^\rho\f$
   */
  vector<double> c12Rho_;
  
  /**
   *  Coupling \f$C_1^\omega+C_2^\omega\f$
   */
  vector<double> c12Omega_;
  
  /**
   *  Coupling \f$C_1^\phi+C_2^\phi\f$
   */
  vector<double> c12Phi_;
  
  /**
   *  Coupling \f$C_2^\rho\f$
   */
  vector<double> c2Rho_;
  
  /**
   *  Coupling \f$C_2^\omega\f$
   */
  vector<double> c2Omega_;
  
  /**
   *  Coupling \f$C_2^\phi\f$
   */
  vector<double> c2Phi_;
  
  /**
   *  Magnetic moments
   */
  vector<double> mu_;
  //@}
};

}

#endif /* HERWIG_KornerKurodaFormFactor_H */
