// -*- C++ -*-
#ifndef HERWIG_LHModel_H
#define HERWIG_LHModel_H
//
// This is the declaration of the LHModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
#include "LHModel.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHModel class is the main class for the implementation of the Little Higgs model in Herwig.
 * In inherits from the StandardModel class and implements the calcuation of the couplings 
 * and masses in the Little Higgs model and storage of the additional couplings.
 *
 * @see \ref LHModelInterfaces "The interfaces"
 * defined for LHModel.
 */
class LHModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  LHModel();

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
  double lambda1() const {return _lambda1;}

  /**
   *  The \f$\lambda_2\f$ top Yukawa coupling
   */
  double lambda2() const {return _lambda2;}

  /**
   *  The sine of the \f$\theta\f$ mixing angle
   */
  double sinTheta() const {return _s;}

  /**
   *  The cosine of the \f$\theta\f$ mixing angle
   */
  double cosTheta() const {return _c;}

  /**
   *  The sine of the \f$\theta'\f$ mixing angle
   */
  double sinThetaPrime() const {return _sp;}

  /**
   *  The cosine of the \f$\theta'\f$ mixing angle
   */
  double cosThetaPrime() const {return _cp;}

  /**
   *  The sine of the Higgs mixing angle
   */
  double sinTheta0() const {return _s0;}

  /**
   *  The cosine of the Higgs mixing angle
   */
  double cosTheta0() const {return _c0;}

  /**
   *  The sine of the pseudoscalar Higgs mixing angle
   */
  double sinThetaP() const {return _sP;}

  /**
   *  The cosine of the pseudoscalar Higgs mixing angle
   */
  double cosThetaP() const {return _cP;}

  /**
   *  The sine of the charged Higgs mixing angle
   */
  double sinThetaPlus() const {return _sPlus;}

  /**
   *  The cosine of the charged Higgs mixing angle
   */
  double cosThetaPlus() const {return _cPlus;}

  /**
   *  The vacuum expection value
   */
  Energy vev() const {return _v;}

  /**
   *  The vacuum expection value
   */
  Energy vevPrime() const {return _v*_vacratio;}

  /**
   *  The \f$f\f$ scale of the non-linear \f$\sigma\f$-model
   */
  Energy f() const {return _f;}
  //@}

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
   * Initialize this object after the setup phase before saving an
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
  LHModel & operator=(const LHModel &);

private:

  /**
   *  Parameters for the model
   */
  //@{
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

  /**
   *  The sine of the pseudoscalar Higgs mixing angle
   */
  double _sP;

  /**
   *  The cosine of the pseudoscalar Higgs mixing angle
   */
  double _cP;

  /**
   *  The sine of the charged Higgs mixing angle
   */
  double _sPlus;

  /**
   *  The cosine of the charged Higgs mixing angle
   */
  double _cPlus;
  //@}

  /**
   *  Additional vertices
   */
  //@{
  /**
   *  WHH Vertex
   */
  AbstractVSSVertexPtr WHHVertex_;
  //@}
};

}

#endif /* HERWIG_LHModel_H */
