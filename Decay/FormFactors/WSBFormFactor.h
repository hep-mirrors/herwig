// -*- C++ -*-
//
// WSBFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WSBFormFactor_H
#define HERWIG_WSBFormFactor_H
//
// This is the declaration of the WSBFormFactor class.
//
#include "ScalarFormFactor.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The WSBFormFactor class is the implementation of the form factor
   *  model of Z.Phys. C29, 637 for the semi-leptonic form factors. It includes
   *  form factors for a number of \f$D\f$, \f$B\f$ and \f$D_s\f$ decays. 
   *  In practice the parameters
   *  of the model were taken from Z.Phys. C34, 103 which includes a number of
   *  decay modes which were not considered in the original paper. 
   *
   *  This form factor model is included both to give an alternative for many modes
   *  to the ISGW models and for use in the factorisation approxmation for hadronic
   *  decays.
   * 
   * @see ScalarFormFactor
   * @see ISGWFormFactor
   * @see ISGW2FormFactor
   *
   */

class WSBFormFactor: public ScalarFormFactor {

public:

  /**
   * Default constructor
   */
  WSBFormFactor();

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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<WSBFormFactor> initWSBFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  WSBFormFactor & operator=(const WSBFormFactor &);

private:

  /** @name Parameters for the form factors */
  //@{
  /**
   * The form factor at \f$q^2=0\f$ for scalar decays.
   */
  vector<double> _F0;

  /**
   * the form factor \f$V\f$ at \f$q^2=0\f$ for vector decays.
   */
  vector<double> _V;

  /**
   * the form factor \f$A_0\f$ at \f$q^2=0\f$ for vector decays.
   */
  vector<double> _A0;

  /**
   * the form factor \f$A_1\f$ at \f$q^2=0\f$ for vector decays.
   */
  vector<double> _A1;

  /**
   * the form factor \f$A_2\f$ at \f$q^2=0\f$ for vector decays.
   */
  vector<double> _A2;

  /**
   * Spin-0 mass for the scalar form factors
   */
  vector<Energy> _mS0;

  /**
   * Spin-1 mass for the scalar form factors
   */
  vector<Energy> _mS1;

  /**
   * Spin-0 mass for the vector form factors
   */
  vector<Energy> _mV0;

  /**
   * Spin-1 mass for the vector form factors
   */
  vector<Energy> _mV1;
  //@}

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * WSBFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::WSBFormFactor,1> {
  /** Typedef of the base class of WSBFormFactor. */
  typedef Herwig::ScalarFormFactor NthBase;
};
/**
 * This template specialization informs ThePEG about the name of the
 * WSBFormFactor class.
 */
template <>
 struct ClassTraits<Herwig::WSBFormFactor>
  : public ClassTraitsBase<Herwig::WSBFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::WSBFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

#endif /* HERWIG_WSBFormFactor_H */
