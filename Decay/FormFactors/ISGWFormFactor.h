// -*- C++ -*-
//
// ISGWFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ISGWFormFactor_H
#define HERWIG_ISGWFormFactor_H
//
// This is the declaration of the <ISGWFormFactor class.
//
#include "ScalarFormFactor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The ISGWFormFactor class is the implementation of the ISGW model of 
   *  Phys. Rev. D39, 799 (1989) for the form-factors. It inherits from the 
   *  ScalarFormFactor class and members
   *  for the calculation of the relevant form factors.
   * 
   * @see ScalarFormFactor
   * @see ISGW2FormFactor
   */

class ISGWFormFactor: public ScalarFormFactor {

public:

  /**
   * Default constructor
   */
  ISGWFormFactor();

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
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,Energy m0,
				      Energy m1,Complex & f0,Complex & fp) const;

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
   * The form factor for the weak decay of a scalar to a tensor.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param h  The form-factor \f$h\f$.
   * @param k  The form-factor \f$k\f$. 
   * @param bp The form-factor \f$b_+\f$.
   * @param bm The form-factor \f$b_-\f$.
   */
  virtual void ScalarTensorFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,
				      Energy m1, complex<InvEnergy2> & h,
				      Complex & k, complex<InvEnergy2> & bp,
				      complex<InvEnergy2> & bm) const;
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

  /** The member which implements all the different form-factors
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f1 The first  form-factor.
   * @param f2 The second form-factor. 
   * @param f3 The third  form-factor.
   * @param f4 The fourth form-factor.
   */
  void formFactor(Energy2 q2,unsigned int iloc,int id0,int id1,Energy m0,
		  Energy m1,Complex & f1,Complex & f2,
		  Complex & f3,Complex & f4) const;
  // general member to calculate all the form-factors

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
  static ClassDescription<ISGWFormFactor> initISGWFormFactor;
  /**
   * Private and non-existent assignment operator.
   */
  ISGWFormFactor & operator=(const ISGWFormFactor &);

private:

  /**
   * The relativistic compensation factor
   */
  double _kappa;

  /** @name Quark masses */
  //@{
  /**
   * The down quark mass
   */
  Energy _mdown;

  /**
   * The up quark mass
   */
  Energy _mup;

  /**
   * The strange quark mass
   */
  Energy _mstrange;

  /**
   * The charm quark mass
   */
  Energy _mcharm;

  /**
   * The bottom quark mass
   */
  Energy _mbottom;

  /**
   * The masses of the quarks as a vector
   */
  vector<Energy> _mquark;
  //@}


  /** @name Wave function parameters */
  //@{
  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{d}\f$ 
   */
  Energy _betaSud;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{s}\f$ 
   */
  Energy _betaSus;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{c}\f$ 
   */
  Energy _betaSuc;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{b}\f$ 
   */
  Energy _betaSub;

  /**
   * The s-wave variational parameters as a vector.
   */
  vector<vector<Energy> > _betaS;

  /**
   * The wavefunction p-wave \f$\beta\f$ variational parameters for  \f$u\bar{d}\f$ 
   */
  Energy _betaPud;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{s}\f$ 
   */
  Energy _betaPus;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{c}\f$ 
   */
  Energy _betaPuc;

  /**
   * The p-wave variational parameters as a vector
   */
  vector<vector<Energy> > _betaP;
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
 * ISGWFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::ISGWFormFactor,1> {
  /** Typedef of the base class of ISGWFormFactor. */
  typedef Herwig::ScalarFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ISGWFormFactor class.
 */
template <>
struct ClassTraits<Herwig::ISGWFormFactor>
  : public ClassTraitsBase<Herwig::ISGWFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::ISGWFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ISGWFormFactor_H */
