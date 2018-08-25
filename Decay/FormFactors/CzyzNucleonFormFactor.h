// -*- C++ -*-
#ifndef Herwig_CzyzNucleonFormFactor_H
#define Herwig_CzyzNucleonFormFactor_H
//
// This is the declaration of the CzyzNucleonFormFactor class.
//

#include "BaryonFormFactor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The CzyzNucleonFormFactor class implements the model of 
 * Phys.Rev. D90 (2014) no.11, 114021 for the proton and neutron form factors
 *
 * @see \ref CzyzNucleonFormFactorInterfaces "The interfaces"
 * defined for CzyzNucleonFormFactor.
 */
class CzyzNucleonFormFactor: public BaryonFormFactor {

public:

  /**
   * The default constructor.
   */
  CzyzNucleonFormFactor();

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
  CzyzNucleonFormFactor & operator=(const CzyzNucleonFormFactor &);

private:

  /**
   *  The masses of the \f$\rho\f$ resonances
   */
  vector<Energy> rhoMasses_;
  
  /**
   *  The widths of the \f$\rho\f$ resonances
   */
  vector<Energy> rhoWidths_;
  
  /**
   *  The masses of the \f$\omega\f$ resonances
   */
  vector<Energy> omegaMasses_;
  
  /**
   *  The widths of the \f$\omega\f$ resonances
   */
  vector<Energy> omegaWidths_;

  /**
   *   The real part of the \f$c_1\f$ couplings
   */
  vector<double> c1Re_;

  /**
   *   The imaginary part of the \f$c_1\f$ couplings
   */
  vector<double> c1Im_;

  /**
   *   The real part of the \f$c_2\f$ couplings
   */
  vector<double> c2Re_;

  /**
   *   The imaginary part of the \f$c_2\f$ couplings
   */
  vector<double> c2Im_;

  /**
   *   The real part of the \f$c_3\f$ couplings
   */
  vector<double> c3Re_;

  /**
   *   The imaginary part of the \f$c_3\f$ couplings
   */
  vector<double> c3Im_;
  
  /**
   *   The real part of the \f$c_4\f$ couplings
   */
  vector<double> c4Re_;

  /**
   *   The imaginary part of the \f$c_4\f$ couplings
   */
  vector<double> c4Im_;

  /**
   *  The complex \f$c_1\f$ coupling
   */
  vector<Complex> c1_;

  /**
   *  The complex \f$c_2\f$ coupling
   */
  vector<Complex> c2_;

  /**
   *  The complex \f$c_3\f$ coupling
   */
  vector<Complex> c3_;

  /**
   *  The complex \f$c_4\f$ coupling
   */
  vector<Complex> c4_;

  /**
   *   Proton magentic moment
   */
  double mup_;

  /**
   *  Neutron magentic moment
   */
  double mun_;
  
  /**
   *   \f$a\f$ parameter
   */
  double a_;

  /**
   *   \f$b\f$ parameter
   */
  double b_;
};

}

#endif /* Herwig_CzyzNucleonFormFactor_H */
