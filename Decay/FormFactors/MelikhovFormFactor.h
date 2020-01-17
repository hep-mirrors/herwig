// -*- C++ -*-
//
// MelikhovFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MelikhovFormFactor_H
#define HERWIG_MelikhovFormFactor_H
//
// This is the declaration of the MelikhovFormFactor class.
//

#include "ScalarFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The MelikhovFormFactor class implements the model of  Phys. Lett. B 380 (1996) 363
 * for the form factors for \f$B\to\pi,\rho\f$. 
 *
 * @see ScalarFormFactor
 *
 */
class MelikhovFormFactor: public ScalarFormFactor {

public:

  /**
   * The default constructor.
   */
  MelikhovFormFactor();

  /** @name Form-Factors */
  //@{
  /**
   * The form factor for the weak decay of a scalar to a scalar. 
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f0 The form-factor \f$f_0\f$. 
   * @param fp The form-factor \f$f_+\f$.
   */
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1,
				      Complex & f0,Complex & fp) const;

  /**
   * The form factor for the weak decay of a scalar to a vector.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param V  The form-factor \f$V\f$
   * @param A0 The form-factor \f$A_0\f$
   * @param A1 The form-factor \f$A_1\f$
   * @param A2 The form-factor \f$A_2\f$
   */
  virtual void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
				      Energy m0, Energy m1, Complex & V,
				      Complex & A0,Complex & A1,Complex & A2) const;
  //@}

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
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MelikhovFormFactor & operator=(const MelikhovFormFactor &) = delete;

private:

  /** @name Parameters for the form factors */
  //@{
  /**
   *  The set of fit parameters to use
   */
  unsigned int _ifit;

  /**
   *  The value of \f$R_+(0)\f$ for the \f$B\to\pi\f$ form factor.
   */
  double _Rplus0;

  /**
   *  The value of \f$M_+\f$ for the \f$B\to\pi\f$ form factor.
   */
  Energy _Mplus;

  /**
   *  The value of \f$n_+\f$ for the \f$B\to\pi\f$ form factor.
   */
  double _nplus;

  /**
   *  The value of \f$R_V(0)\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _RV0;

  /**
   *  The value of \f$M_V\f$ for the \f$B\to\rho\f$ form factor.
   */
  Energy _MV;

  /**
   *  The value of \f$n_V\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _nV;

  /**
   *  The value of \f$R_1(0)\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _R10;

  /**
   *  The value of \f$M_1\f$ for the \f$B\to\rho\f$ form factor.
   */
  Energy _M1;

  /**
   *  The value of \f$n_1\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _n1;

  /**
   *  The value of \f$R_2(0)\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _R20;

  /**
   *  The value of \f$M_2\f$ for the \f$B\to\rho\f$ form factor.
   */
  Energy _M2;

  /**
   *  The value of \f$n_2\f$ for the \f$B\to\rho\f$ form factor.
   */
  double _n2;
  //@}
};

}

#endif /* HERWIG_MelikhovFormFactor_H */
