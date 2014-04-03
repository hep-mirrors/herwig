// -*- C++ -*-
//
// DipoleIOperator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleIOperator_H
#define HERWIG_DipoleIOperator_H
//
// This is the declaration of the DipoleIOperator class.
//

#include "Herwig++/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief DipoleIOperator implements the I(\epsilon)
 * insertion operator.
 *
 */
class DipoleIOperator: public MatchboxInsertionOperator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleIOperator();

  /**
   * The destructor.
   */
  virtual ~DipoleIOperator();
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
   * Vector to contain heavy flavour id's
   * n_F = NHeavy.size()
   */
  vector<int> NHeavy;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleIOperator & operator=(const DipoleIOperator &);

};

}

#endif /* HERWIG_DipoleIOperator_H */
