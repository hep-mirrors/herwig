// -*- C++ -*-
#ifndef THEPEG_LHTPModel_H
#define THEPEG_LHTPModel_H
//
// This is the declaration of the LHTPModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "LHTPModel.fh"

namespace Herwig {

/**
 * The LHTPModel class is the main class for the
 * implementation of the Little Higgs model with T-parity
 *
 * @see \ref LHTPModelInterfaces "The interfaces"
 * defined for LHTPModel.
 */
class LHTPModel: public BSMModel {

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
  Energy vev() const { return v_; }

  /**
   *  The \f$f\f$ scale of the non-linear \f$\sigma\f$-model
   */
  Energy f() const { return f_; }

  /**
   *  \f$\sin\alpha\f$
   */
  double sinAlpha() const { return salpha_; }

  /**
   *  \f$\cos\alpha\f$
   */
  double cosAlpha() const { return calpha_; }

  /**
   *  \f$\sin\beta\f$
   */
  double sinBeta() const { return sbeta_; }

  /**
   *  \f$\cos\beta\f$
   */
  double cosBeta() const { return cbeta_; }

  /**
   *  \f$\sin\theta_H\f$
   */
  double sinThetaH() const { return sthetaH_; }

  /**
   *  \f$\cos\theta_H\f$
   */
  double cosThetaH() const { return cthetaH_; }

  /**
   *  \f$\sin\theta_L\f$
   */
  double sinThetaL() const { return sL_;}

  /**
   *  \f$\cos\theta_L\f$
   */
  double cosThetaL() const { return cL_;}

  /**
   *  \f$\sin\theta_R\f$
   */
  double sinThetaR() const { return sR_;}

  /**
   *  \f$\cos\theta_R\f$
   */
  double cosThetaR() const { return cR_;}
  //@}

  /**
   *  Yukwawa for T-odd fermions
   */
  //@{
  /**
   *  The \f$\kappa_q\f$ parameter which controls the properties of the
   *  T-odd quarks
   */
  double kappaQuark() const { return kappaQuark_;}

  /**
   *  The \f$\kappa_\ell\f$ parameter which controls the properties of the
   *  T-odd leptons
   */
  double kappaLepton() const { return kappaLepton_;}
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
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPModel & operator=(const LHTPModel &) = delete;

private:

  /**
   *  The constant for the non-linear \f$\sigma\f$ model
   */
  Energy f_;

  /**
   *  @name The mixing in the top quark sector
   */
  //@{
  /**
   *  \f$\sin\alpha\f$, taken as an input
   */
  double salpha_;

  /**
   *  \f$\cos\alpha\f$
   */
  double calpha_;

  /**
   *  \f$\sin\beta\f$
   */
  double sbeta_;

  /**
   *  \f$\cos\beta\f$
   */
  double cbeta_;
  //@}

  /**
   *  @name Mixing of the heavy photon and Z
   */
  //@{
  /**
   *  \f$\sin\theta_H\f$
   */
  double sthetaH_;

  /**
   *  \f$\cos\theta_H\f$
   */
  double cthetaH_;
  //@}

  /**
   *  @name Mixings in the top sector
   */
  //@{
  /**
   *  \f$\sin\theta_L\f$
   */
  double sL_;

  /**
   *  \f$\cos\theta_L\f$
   */
  double cL_;

  /**
   *  \f$\sin\theta_R\f$
   */
  double sR_;

  /**
   *  \f$\cos\theta_R\f$
   */
  double cR_;
  //@}

  /**
   *  The \f$\kappa_q\f$ parameter which controls the properties of the
   *  T-odd quarks
   */
  double kappaQuark_;

  /**
   *  The \f$\kappa_\ell\f$ parameter which controls the properties of the
   *  T-odd leptons
   */
  double kappaLepton_;

  /**
   *  The mass of the Standard Model higgs
   */
  Energy mh_;

  /**
   *  The vacuum expection valve
   */
  Energy v_;

  /**
   *  The \f$g\f$ coupling
   */
  double g_;

  /**
   *  the \f$g'\f$ coupling
   */
  double gp_;

  /**
   *  Method for evaluating the masses
   */
  bool approximate_;

  /**
   *  WHH Vertex
   */
  AbstractVSSVertexPtr WHHVertex_;
};

}

#endif /* THEPEG_LHTPModel_H */
