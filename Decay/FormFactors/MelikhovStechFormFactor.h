// -*- C++ -*-
//
// MelikhovStechFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MelikhovStechFormFactor_H
#define HERWIG_MelikhovStechFormFactor_H
//
// This is the declaration of the MelikhovStechFormFactor class.
//

#include "ScalarFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The MelikhovStechFormFactor class is the implementation of the form factors
 * from Phys. Rev. D62  014006 (2000).
 *
 * @see ScalarFormFactor
 */
class MelikhovStechFormFactor: public ScalarFormFactor {

public:

  /**
   * The default constructor.
   */
  MelikhovStechFormFactor();

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

  /**
   * The form factor for the weak penguin decay of a scalar meson to a scalar meson.
   * This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * scalar penguin form factors. 
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param fT The form factor \f$f_T\f$.
   */
  virtual void ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
					   Energy m0, Energy m1,Complex & fT) const;

  /**
   * The form factor for the weak penguin decay of a scalar meson to a vector meson. 
   * This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * vector penguin form factors.
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
  virtual void ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MelikhovStechFormFactor> initMelikhovStechFormFactor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MelikhovStechFormFactor & operator=(const MelikhovStechFormFactor &) = delete;

private:


  /** @name Parameters for the scalar to scalar form factors */
  //@{

  /**
   *  The value of \f$F_+(0)\f$ for the form factors.
   */
  vector<double> _fplus0;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$F_+\f$ form factor.
   */
  vector<double> _sigma1fp;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$F_+\f$ form factor.
   */
  vector<double> _sigma2fp;

  /**
   *  The value of \f$F_0(0)\f$ for the form factors.
   */
  vector<double> _f00;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$F_0\f$ form factor.
   */
  vector<double> _sigma1f0;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$F_0\f$ form factor.
   */
  vector<double> _sigma2f0;

  /**
   *  The value of \f$F_T(0)\f$ for the form factors.
   */
  vector<double> _fT0;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$F_T\f$ form factor.
   */
  vector<double> _sigma1fT;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$F_T\f$ form factor.
   */
  vector<double> _sigma2fT;

  //@}

  /** @name Parameters for the scalar to vector form factors */
  //@{
  /**
   *  The value of \f$V(0)\f$ for the form factors
   */
  vector<double> _V0;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$V_0\f$ form factor.
   */
  vector<double> _sigma1V0;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$V_0\f$ form factor.
   */
  vector<double> _sigma2V0;

  /**
   *  The value of \f$A_0(0)\f$ for the form factors
   */
  vector<double> _A00;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$A_0\f$ form factor.
   */
  vector<double> _sigma1A0;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$A_0\f$ form factor.
   */
  vector<double> _sigma2A0;

  /**
   *  The value of \f$A_1(0)\f$ for the form factors
   */
  vector<double> _A10;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$A_1\f$ form factor.
   */
  vector<double> _sigma1A1;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$A_1\f$ form factor.
   */
  vector<double> _sigma2A1;

  /**
   *  The value of \f$A_2(0)\f$ for the form factors
   */
  vector<double> _A20;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$A_2\f$ form factor.
   */
  vector<double> _sigma1A2;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$A_2\f$ form factor.
   */
  vector<double> _sigma2A2;

  /**
   *  The value of \f$T_1(0)\f$ for the form factors
   */
  vector<double> _T10;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$T_1\f$ form factor.
   */
  vector<double> _sigma1T1;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$T_1\f$ form factor.
   */
  vector<double> _sigma2T1;

  /**
   *  The value of \f$T_2(0)\f$ for the form factors
   */
  vector<double> _T20;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$T_2\f$ form factor.
   */
  vector<double> _sigma1T2;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$T_2\f$ form factor.
   */
  vector<double> _sigma2T2;

  /**
   *  The value of \f$T_3(0)\f$ for the form factors
   */
  vector<double> _T30;

  /**
   *  The \f$\sigma_1\f$ parameter for the \f$T_2\f$ form factor.
   */
  vector<double> _sigma1T3;

  /**
   *  The \f$\sigma_2\f$ parameter for the \f$T_2\f$ form factor.
   */
  vector<double> _sigma2T3;
  //@}

  /**
   *   The scalar mass for the \f$q^2\f$ dependence of the form factor.
   */
  vector<Energy> _massP;

  /**
   *   The vector mass for the \f$q^2\f$ dependence of the form factor.
   */
  vector<Energy> _massV;

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MelikhovStechFormFactor. */
template <>
 struct BaseClassTrait<Herwig::MelikhovStechFormFactor,1> {
  /** Typedef of the first base class of MelikhovStechFormFactor. */
  typedef Herwig::ScalarFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MelikhovStechFormFactor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MelikhovStechFormFactor>
  : public ClassTraitsBase<Herwig::MelikhovStechFormFactor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MelikhovStechFormFactor"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MelikhovStechFormFactor class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MelikhovStechFormFactor_H */
