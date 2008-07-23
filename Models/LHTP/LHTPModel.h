// -*- C++ -*-
#ifndef THEPEG_LHTPModel_H
#define THEPEG_LHTPModel_H
//
// This is the declaration of the LHTPModel class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "LHTPModel.fh"

namespace Herwig {

/**
 * The LHTPModel class is the main class for the
 * implementation of the Little Higgs model with T-parity
 *
 * @see \ref LHTPModelInterfaces "The interfaces"
 * defined for LHTPModel.
 */
class LHTPModel: public StandardModel {

public:

  /**
   * The default constructor.
   */
  LHTPModel();

  /**
   *  Access to the parameters of the model
   */
  //@{
  /**
   *  The vacuum expection value
   */
  Energy vev() const { return _v; }

  /**
   *  The \f$f\f$ scale of the non-linear \f$\sigma\f$-model
   */
  Energy f() const { return _f; }

  /**
   *  \f$\sin\alpha\f$
   */
  double sinAlpha() const { return _salpha; }

  /**
   *  \f$\cos\alpha\f$
   */
  double cosAlpha() const { return _calpha; }

  /**
   *  \f$\sin\beta\f$
   */
  double sinBeta() const { return _sbeta; }

  /**
   *  \f$\cos\beta\f$
   */
  double cosBeta() const { return _cbeta; }

  /**
   *  \f$\sin\theta_H\f$
   */
  double sinThetaH() const { return _sthetaH; }

  /**
   *  \f$\cos\theta_H\f$
   */
  double cosThetaH() const { return _cthetaH; }
  //@}

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
   *  Reset the mass of a ParticleData object
   */
  void resetMass(long id, Energy mass);

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
   *  Calculate the mixing in the top sector of the model
   *  and the masses of the T-odd and T-even heavy tops
   *  The mixings are calculated by solving Eqns 2.22 and 2.24 of hep-ph/0506042
   *  for \f$\lambda_1$ and \f$\lambda_2\f$ given the input value of \f$\sin\alpha\f$
   *  and the top mass.
   */
  void topMixing(Energy & MTp, Energy & MTm);

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LHTPModel> initLHTPModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPModel & operator=(const LHTPModel &);

private:

  /**
   *  The constant for the non-linear \f$\sigma\f$ model
   */
  Energy _f;

  /**
   *  @name The mixing in the top quark sector
   */
  //@{
  /**
   *  \f$\sin\alpha\f$, taken as an input
   */
  double _salpha;

  /**
   *  \f$\cos\alpha\f$
   */
  double _calpha;

  /**
   *  \f$\sin\beta\f$
   */
  double _sbeta;

  /**
   *  \f$\cos\beta\f$
   */
  double _cbeta;
  //@}

  /**
   *  @name Mixing of the heavy photon and Z
   */
  //@{
  /**
   *  \f$\sin\theta_H\f$
   */
  double _sthetaH;

  /**
   *  \f$\cos\theta_H\f$
   */
  double _cthetaH;
  //@}

  /**
   *  The \f$\kappa\f$ parameter which controls the properties of the
   *  T-odd fermions
   */
  double _kappa;

  /**
   *  The mass of the Standard Model higgs
   */
  Energy _mh;

  /**
   *  The vacuum expection valve
   */
  Energy _v;

  /**
   *  The \f$g\f$ coupling
   */
  double _g;

  /**
   *  the \f$g'\f$ coupling
   */
  double _gp;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHTPModel. */
template <>
struct BaseClassTrait<Herwig::LHTPModel,1> {
  /** Typedef of the first base class of LHTPModel. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHTPModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHTPModel>
  : public ClassTraitsBase<Herwig::LHTPModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHTPModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHTPModel is implemented. It may also include several, space-separated,
   * libraries if the class LHTPModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHTPModel.so"; }
};

/** @endcond */

}

#endif /* THEPEG_LHTPModel_H */
