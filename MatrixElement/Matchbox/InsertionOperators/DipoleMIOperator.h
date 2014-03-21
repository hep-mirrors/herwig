// -*- C++ -*-
//
// DipoleMIOperator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleMIOperator_H
#define HERWIG_DipoleMIOperator_H
//
// This is the declaration of the DipoleMIOperator class.
//

#include "Herwig++/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Martin Stoll, Christian Reuschle
 *
 * \brief DipoleMIOperator implements the massive I operator.
 *
 * It is lim_{m->0}I(m)->I(m=0)
 *
 * More precisely: DipoleIOperator applies only for massless
 * partons while DipoleMIOperator applies for massless AND
 * massive partons. In the sum over emitter and spectator,
 * if both are massless, the difference between the massive
 * and massless I operator is returned. However, if at least
 * one is massive, the massive I operator is returned. If both
 * are incoming, then there will be no associated contribution
 * from DipoleMIOperator.
 *
 * In the case of g->QQbar, which is emitter=g and maybe mass-
 * less spectator, the additional terms through the sums over
 * massive flavours should be correctly accounted for.
 *
 * DipoleMIOperator does only apply in the expanded convention,
 * and also not for dimensional reduction.
 *
 * DipoleMIOperator trusts that initial state particles are not
 * massive.
 *
 */
class DipoleMIOperator: public MatchboxInsertionOperator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleMIOperator();

  /**
   * The destructor.
   */
  virtual ~DipoleMIOperator();
  //@}

public:

  /**
   * Set the XComb object steering the Born matrix
   * element this class represents virtual corrections to.
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual double me2() const;

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const;

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const;

//   /**
//    * Return true, if contributions exist to
//    * the given parton pair.
//    */
//   bool apply(tcPDPtr, tcPDPtr) const;

  /**
   * Return true, if contributions exist to
   * the given parton.
   */
  bool apply(tcPDPtr) const;

  /**
   * Return true, if this virtual correction
   * applies to the given process.
   */
  virtual bool apply(const cPDVector&) const;

public:
  
  /**
   * Triangular / Kallen function
   */
  template <class T>
  inline T rootOfKallen (T a, T b, T c) const {
    return sqrt( a*a + b*b + c*c - 2.*( a*b+a*c+b*c ) ); }

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


private:

  /**
   * C_A
   */
  double CA;

  /**
   * C_F
   */
  double CF;

  /**
   * \gamma_q
   */
  double gammaQuark;

  /**
   * \gamma_g
   */
  double gammaGluon;
  
  /**
   * \beta_0
   * The Matchbox convention is \beta_0=\gamma_g with the counterterm part as n/\epsilon*\beta_0.
   * However, usually \beta_0 is defined as \beta_0=2*\gamma_g, and instead the counterterm part as n/\epsilon*1/2*\beta_0.
   */
  double betaZero;

  /**
   * \Gamma_q, finite term
   */
  double GammaQuark(const ParticleData&,Energy2) const;
  
  /**
   * \Gamma_g, finite term
   */
  double GammaGluon() const;

  /**
   * K_q
   */
  double KQuark;

  /**
   * K_g
   */
  double KGluon;
  
  /**
   * V_j
   */
  double Vj(const ParticleData&, const ParticleData&, Energy2,double,bool=false) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleMIOperator & operator=(const DipoleMIOperator &);

};

}

#endif /* HERWIG_DipoleMIOperator_H */
