// -*- C++ -*-
//
// BallZwickyVectorFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BallZwickyVectorFormFactor_H
#define HERWIG_BallZwickyVectorFormFactor_H
//
// This is the declaration of the BallZwickyVectorFormFactor class.
//
#include "ScalarFormFactor.h"
namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The BallZwickyVectorFormFactor class implements the form-factors
   *  of hep-ph/0412079 for the B meson to light vector mesons.
   *
   *  This class is one of the few which includes the penguin form factors in addition
   *  to the standard weak decay form factors.
   *
   * @see ScalarFormFactor
   * @see BallZwickyScalarFormFactor
   */

class BallZwickyVectorFormFactor: public ScalarFormFactor {

public:

  /**
   * Default constructor
   */
  BallZwickyVectorFormFactor();

  /** @name Form-Factors */
  //@{
  /**
   * The form factor for the weak decay of a scalar to a vector.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param A0 The form factor \f$A_0\f$
   * @param A1 The form factor \f$A_1\f$
   * @param A2 The form factor \f$A_2\f$
   * @param V  The form factor \f$V\f$
   */
  void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
			      Energy m0, Energy m1,Complex & A0,
			      Complex & A1,Complex & A2, Complex & V) const;

  /**
   * The form factor for the weak penguin decay of a scalar meson to a vector meson. 
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param T1 The form factor \f$T_1\f$.
   * @param T2 The form factor \f$T_2\f$.
   * @param T3 The form factor \f$T_3\f$.
   */
  void ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				   Energy m0, Energy m1, Complex & T1,
				   Complex & T2, Complex & T3) const;
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
  BallZwickyVectorFormFactor & operator=(const BallZwickyVectorFormFactor &) = delete;

private:

  /** @name Coefficients for the form factors.*/
  //@{

  /**
   * The coefficient \f$r_1\f$ for the \f$V(q^2)\f$ form factor.
   */
  vector<double> _Vr1;

  /**
   * The coefficient \f$r_2\f$ for the \f$V(q^2)\f$ form factor.
   */
  vector<double> _Vr2;

  /**
   * The coefficient \f$r_1\f$ for the \f$A_0(q^2)\f$ form factor.
   */
  vector<double> _A0r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$A_0(q^2)\f$ form factor.
   */
  vector<double> _A0r2;

  /**
   * The coefficient \f$r_1\f$ for the \f$A_1(q^2)\f$ form factor.
   */
  vector<double> _A1r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$A_1(q^2)\f$ form factor.
   */
  vector<double> _A1r2;

  /**
   * The coefficient \f$r_1\f$ for the \f$A_2(q^2)\f$ form factor.
   */
  vector<double> _A2r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$A_2(q^2)\f$ form factor.
   */
  vector<double> _A2r2;

  /**
   * The coefficient \f$r_1\f$ for the \f$T_1(q^2)\f$ form factor.
   */
  vector<double> _T1r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$T_1(q^2)\f$ form factor.
   */
  vector<double> _T1r2;

  /**
   * The coefficient \f$r_1\f$ for the \f$T_2(q^2)\f$ form factor.
   */
  vector<double> _T2r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$T_2(q^2)\f$ form factor.
   */
  vector<double> _T2r2;

  /**
   * The coefficient \f$r_1\f$ for the \f$\tilde{T}_3(q^2)\f$ form factor.
   */
  vector<double> _T3r1;

  /**
   * The coefficient \f$r_2\f$ for the \f$\tilde{T}_3(q^2)\f$ form factor.
   */
  vector<double> _T3r2;
  // the constants for the form-factors
  //@}

  /** @name Masses for the form factors.*/
  //@{

  /**
   * The mass \f$m_R^2\f$ for the \f$V(q^2)\f$ form factor.
   */
  vector<Energy2> _VmR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$V(q^2)\f$ form factor.
   */
  vector<Energy2> _Vmfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$A_0(q^2)\f$ form factor.
   */
  vector<Energy2> _A0mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$A_0(q^2)\f$ form factor.
   */
  vector<Energy2> _A0mfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$A_1(q^2)\f$ form factor.
   */
  vector<Energy2> _A1mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$A_1(q^2)\f$ form factor.
   */
  vector<Energy2> _A1mfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$A_2(q^2)\f$ form factor.
   */
  vector<Energy2> _A2mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$A_2(q^2)\f$ form factor.
   */
  vector<Energy2> _A2mfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$T_1(q^2)\f$ form factor.
   */
  vector<Energy2> _T1mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$T_1(q^2)\f$ form factor.
   */
  vector<Energy2> _T1mfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$T_2(q^2)\f$ form factor.
   */
  vector<Energy2> _T2mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$T_2(q^2)\f$ form factor.
   */
  vector<Energy2> _T2mfit2;

  /**
   * The mass \f$m_R^2\f$ for the \f$\tilde{T}_3(q^2)\f$ form factor.
   */
  vector<Energy2> _T3mR2;

  /**
   * The mass \f$m_{\rm fit}^2\f$ for the \f$\tilde{T}_3(q^2)\f$ form factor.
   */
  vector<Energy2> _T3mfit2;
  // the masses for the form-factors
  //@}

  /**
   * Cut-off parameter for the switch to a small \f$q^2\f$ expansion for the \f$T_3\f$
   * form factor.
   */
  Energy2 _cutoff;
};

}

#endif /* HERWIG_BallZwickyVectorFormFactor_H */
