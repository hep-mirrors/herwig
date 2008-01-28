// -*- C++ -*-
#ifndef HERWIG_LHModel_H
#define HERWIG_LHModel_H
//
// This is the declaration of the LHModel class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "LittleHiggsModel.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LHModel class.
 *
 * @see \ref LHModelInterfaces "The interfaces"
 * defined for LHModel.
 */
class LHModel: public StandardModel {

public:

  /**
   * The default constructor.
   */
  inline LHModel();

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

public:

  /**
   *  Access to the parameters of the model
   */
  //@{
  /**
   *  The \f$\lambda_1\f$ top Yukawa coupling
   */
  inline double lambda1() const;

  /**
   *  The \f$\lambda_2\f$ top Yukawa coupling
   */
  inline double lambda2() const;

  /**
   *  The sine of the \f$\theta\f$ mixing angle
   */
  inline double sinTheta() const;

  /**
   *  The cosine of the \f$\theta\f$ mixing angle
   */
  inline double cosTheta() const;

  /**
   *  The sine of the \f$\theta'\f$ mixing angle
   */
  inline double sinThetaPrime() const;

  /**
   *  The cosine of the \f$\theta'\f$ mixing angle
   */
  inline double cosThetaPrime() const;

  /**
   *  The sine of the Higgs mixing angle
   */
  inline double sinTheta0() const;

  /**
   *  The cosine of the Higgs mixing angle
   */
  inline double cosTheta0() const;

  /**
   *  The vacuum expection value
   */
  inline Energy vev() const;

  /**
   *  The vacuum expection value
   */
  inline Energy vevPrime() const;

  /**
   *  The \f$f\f$ scale of the non-linear \f$\sigma\f$-model
   */
  inline Energy f() const;
  //@}

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

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
  static ClassDescription<LHModel> initLHModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHModel & operator=(const LHModel &);

private:

  /**
   *  Parameters for the model
   */
  //@{
  /**
   *  The \f$g\f$ coupling
   */
  double _g;

  /**
   *  the \f$g'\f$ coupling
   */
  double _gp;

  /**
   *  The value of \f$\cot\theta\f$ for the mixing with the \f$g\f$ coupling 
   */
  double _cott;

  /**
   *  The value of \f$\tan\theta'\f$ for the mixing with the \f$g'\f$ coupling
   */
  double _tantp;

  /**
   *  The vacuum expection valve
   */
  Energy _v;

  /**
   *  The ratio \f$\lambda_2/\lambda_1\f$
   */
  double _lamratio;

  /**
   *  The mass of the lightest Higgs boson
   */
  Energy _mH;

  /**
   *  The ratio of the vacuum exception values \f$v'/v\f$
   */
  double _vacratio;

  /**
   *  The scale for the non-linear \f$\sigma\f$-model
   */
  Energy _f;

  /**
   *  The top Yukawa couplings
   */
  double _lambda1,_lambda2;

  /**
   *  The sine of the \f$\theta\f$ mixing angle
   */
  double _s;

  /**
   *  The cosine of the \f$\theta\f$ mixing angle
   */
  double _c;

  /**
   *  The sine of the \f$\theta'\f$ mixing angle
   */
  double _sp;

  /**
   *  The cosine of the \f$\theta'\f$ mixing angle
   */
  double _cp;

  /**
   *  The sine of the Higgs mixing angle
   */
  double _s0;

  /**
   *  The cosine of the Higgs mixing angle
   */
  double _c0;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHModel. */
template <>
struct BaseClassTrait<Herwig::LHModel,1> {
  /** Typedef of the first base class of LHModel. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHModel>
  : public ClassTraitsBase<Herwig::LHModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHModel is implemented. It may also include several, space-separated,
   * libraries if the class LHModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#include "LittleHiggsModel.icc"

#endif /* HERWIG_LHModel_H */
