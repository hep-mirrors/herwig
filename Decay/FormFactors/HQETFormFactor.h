// -*- C++ -*-
//
// HQETFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HQETFormFactor_H
#define HERWIG_HQETFormFactor_H
//
// This is the declaration of the HQETFormFactor class.
//

#include "ScalarFormFactor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The HQETFormFactor class implements the form factors for the decay of
 * B mesons into D and \f$D^*\f$ mesons using form factors based on Heavy Quark
 * Effective theory of the parameterisation of hep-ph/9712417 away from the 
 * heavy quark limit.
 *
 * @see \ref HQETFormFactorInterfaces "The interfaces"
 * defined for HQETFormFactor.
 */
class HQETFormFactor: public ScalarFormFactor {

public:

  /**
   * The default constructor.
   */
  HQETFormFactor();

  /** @name Form Factors */
  //@{
  /**
   * The form factor for the weak decay of a scalar to a scalar.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f0 The form factor \f$f_0\f$. 
   * @param fp The form factor \f$f_+\f$.
   */
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1,Complex & f0,
				      Complex & fp) const;

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
  virtual void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
				      Energy m0, Energy m1,Complex & A0,
				      Complex & A1,Complex & A2, Complex & V) const;
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HQETFormFactor> initHQETFormFactor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HQETFormFactor & operator=(const HQETFormFactor &);

private:

  /**
   *  The normalisation parameter for the scalar form factors
   */
  double _f1scalar;

  /**
   *  The normalisation parameter for the vector form factors
   */
  double _f1vector;

  /**
   *  The form factor \f$R_1(\omega)\f$ at \f$\omega=1\f$.
   */
  double _r1;

  /**
   *  The form factor \f$R_2(\omega)\f$ at \f$\omega=1\f$.
   */
  double _r2;

  /**
   *  The slope parameter \f$\rho^2\f$ for the scalar form factors
   */
  double _rho2scalar;

  /**
   *  The slope parameter \f$\rho^2\f$ for the vector form factors
   */
  double _rho2vector;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HQETFormFactor. */
template <>
struct BaseClassTrait<Herwig::HQETFormFactor,1> {
  /** Typedef of the first base class of HQETFormFactor. */
  typedef Herwig::ScalarFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HQETFormFactor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HQETFormFactor>
  : public ClassTraitsBase<Herwig::HQETFormFactor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HQETFormFactor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HQETFormFactor is implemented. It may also include several, space-separated,
   * libraries if the class HQETFormFactor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HQETFormFactor_H */
