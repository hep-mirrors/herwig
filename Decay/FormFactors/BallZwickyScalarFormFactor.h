// -*- C++ -*-
//
// BallZwickyScalarFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BallZwickyScalarFormFactor_H
#define HERWIG_BallZwickyScalarFormFactor_H
//
// This is the declaration of the BallZwickyScalarFormFactor class.
//
#include "ScalarFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This class is the implementation of the form-factors of PRD71 014015 (2005) for
 *  the form-factor for the decay of a B-meson to a light pseudoscalar meson.
 *
 *  This class is one of the few which includes the penguin form factors in addition
 *  to the standard weak decay form factors.
 *
 * @see ScalarFormFactor
 * @see BallZwickyVectorFormFactor
 */
class BallZwickyScalarFormFactor: public ScalarFormFactor {
  
public:

  /**
   * Default constructor
   */
  BallZwickyScalarFormFactor();

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
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,
				      int id1, Energy m0, Energy m1,
				      Complex & f0,Complex & fp) const;

  /**
   * The form factor for the weak penguin decay of a scalar meson to a scalar meson.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param fT The form factor \f$f_T\f$.
   */
  void ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				   Energy m0, Energy m1,Complex & fT) const;
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
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  BallZwickyScalarFormFactor & operator=(const BallZwickyScalarFormFactor &) = delete;

private:

  /** @name Coefficients for the form factors.*/
  //@{

  /**
   * The coefficient \f$r_1\f$ for the \f$f_0(q^2)\f$ form factor.
   */
  vector<double> _r10;

  /**
   * The coefficient \f$r_2\f$ for the \f$f_0(q^2)\f$ form factor.
   */
  vector<double> _r20;

  /**
   * The coefficient \f$r_1\f$ for the \f$f_+(q^2)\f$ form factor.
   */
  vector<double> _r1plus;

  /**
   * The coefficient \f$r_2\f$ for the \f$f_+(q^2)\f$ form factor.
   */
  vector<double> _r2plus;

  /**
   * The coefficient \f$r_1\f$ for the \f$f_T(q^2)\f$ form factor.
   */
  vector<double> _r1T;

  /**
   * The coefficient \f$r_2\f$ for the \f$f_T(q^2)\f$ form factor.
   */
  vector<double> _r2T;
  //@}

  /** @name Masses for the form factors.*/
  //@{

  /**
   * The mass \f$(m_1)^2\f$ \f$f_0(q^2)\f$ form factor.
   */
  vector<Energy2> _m120;

  /**
   * The mass \f$m_{\rm fit}^2\f$ \f$f_0(q^2)\f$ form factor.
   */
  vector<Energy2> _mfit20;

  /**
   * The mass \f$(m_1)^2\f$ \f$f_+(q^2)\f$ form factor.
   */
  vector<Energy2> _m12plus;

  /**
   * The mass \f$m_{\rm fit}^2\f$ \f$f_+(q^2)\f$ form factor.
   */
  vector<Energy2> _mfit2plus;

  /**
   * The mass \f$(m_1)^2\f$ \f$f_T(q^2)\f$ form factor.
   */
  vector<Energy2> _m12T;

  /**
   * The mass \f$m_{\rm fit}^2\f$ \f$f_T(q^2)\f$ form factor.
   */
  vector<Energy2> _mfit2T;
  //@}

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;

};

}

#endif /* HERWIG_BallZwickyScalarFormFactor_H */
