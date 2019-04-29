// -*- C++ -*-
//
// DipoleMIOperator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleMIOperator_H
#define HERWIG_DipoleMIOperator_H
//
// This is the declaration of the DipoleMIOperator class.
//

#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Daniel Rauch, Christian Reuschle,
 *         Martin Stoll
 *
 * \brief DipoleMIOperator implements the I(\epsilon) 
 * insertion operator for the massive case.
 * DipoleMIOperator does only apply for expanded con-
 * vention and also not for dimensional reduction.
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
   * Return true, if this virtual correction
   * applies to the given process.
   */
  virtual bool apply(const cPDVector&) const;

  /**
   * Return true, if contributions exist to
   * the given parton.
   */
  bool apply(tcPDPtr) const;

  /**
   * Return a vector of PDG codes of the light flavours,
   * which are contained in the jet particle group.
   */
  vector<int> NLightJetVec() const;

  /**
   * Return a vector of PDG codes of the heavy flavours,
   * which are contained in the jet particle group.
   */
  vector<int> NHeavyJetVec() const;

  /**
   * Return a vector of PDG codes of the light flavours,
   * which are contained in the associated Born sub-process.
   */
  vector<int> NLightBornVec() const;

  /**
   * Return a vector of PDG codes of the heavy flavours,
   * which are contained in the associated Born sub-process.
   */
  vector<int> NHeavyBornVec() const;

  /**
   * Return a vector of PDG codes of the light flavours,
   * which are contained in the proton particle group.
   */
  vector<int> NLightProtonVec() const;

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
   * The Matchbox convention is \beta_0=\gamma_g. Often however,
   * \beta_0 is defined in the literature as \beta_0=2*\gamma_g.
   * Be aware of consistent usage!
   * In the massive case (see hep-ph/0011222v3):
   *      \beta_0 = 11/3*C_A - 4/3*T_R*(N_f+N_F)
   * with T_R=1/2, N_f the number of light flavours ,and N_F the
   * number of heavy flavours, which originate in the splittings
   * g->qqbar or g->QQbar.
   * In our conventions, however:
   *      \beta_0 = 11/6*C_A - 2/3*T_R*(N_f+N_F)
   * The "massive" \beta_0 applies as soon as we define massive
   * flavours in the jet particle group. Be aware that some OLP
   * might generically(!) use something like N_f=5 and N_F=1 in
   * their definition of \beta_0!
   */
  double betaZero;

  /**
   * K_q
   */
  double KQuark;

  /**
   * K_g
   */
  double KGluon;
  
private:

  /**
   * V_j, non-singular terms
   */
  // double Vj(const ParticleData&, const ParticleData&, Energy2, double, bool=false) const;
  double Vj(const ParticleData&, const ParticleData&, Energy2, double, double, int, bool=false) const;

  /**
   * V^{(s)}, double pole part in expanded convention.
   */
  double VsDoublePole(const ParticleData&, const ParticleData&) const;

  /**
   * V^{(s)}, single pole part in expanded convention.
   */
  double VsSinglePole(const ParticleData&, const ParticleData&, Energy2) const;

  /**
   * \Gamma_q, finite term
   */
  double GammaQuark(const ParticleData&) const;
  
  /**
   * \Gamma_g, finite term
   */
  double GammaGluon() const;

  /**
   * \Gamma_q, single pole term
   */
  double GammaQuarkSinglePole(const ParticleData&) const;
  
  /**
   * \Gamma_g, single pole term
   */
  double GammaGluonSinglePole() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleMIOperator & operator=(const DipoleMIOperator &) = delete;

};

}

#endif /* HERWIG_DipoleMIOperator_H */
