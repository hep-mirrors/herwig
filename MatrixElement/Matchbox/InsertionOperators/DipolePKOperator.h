// -*- C++ -*-
//
// DipolePKOperator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipolePKOperator_H
#define HERWIG_DipolePKOperator_H
//
// This is the declaration of the DipolePKOperator class.
//

#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "ThePEG/PDF/PDF.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Christian Reuschle
 *
 * \brief DipolePKOperator implements the P+K
 * insertion operator.
 *
 */
class DipolePKOperator: public MatchboxInsertionOperator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipolePKOperator();

  /**
   * The destructor.
   */
  virtual ~DipolePKOperator();
  //@}

public:

  /**
   * Set the XComb object steering the Born matrix
   * element this class represents virtual corrections to.
   */
  virtual void setXComb(tStdXCombPtr xc);
  
  /**
   * Set parameters for new alpha parameter.
   */
  virtual void setAlpha (double alpha)const;

  /**
   * Return the number of additional random variables
   * needed to calculate this virtual correction.
   * We treat all integrations on equal footing.
   */
  virtual int nDimAdditional() const { return 1; }

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
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual CrossSection dSigHatDR() const { 
    return sqr(hbarc) * me2() * lastBorn()->lastXComb().jacobian() / (2.*lastSHat());
  }

public:

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

  /**
   * The initial-final contribution
   */
  double gammaSoft() const;

  /**
   * [(1/(1-z))*ln((1-z)/z)]_+
   */
  double softLogByz(tcPDPtr p) const;

  /**
   * [(1/(1-z))*ln((1-z))]_+
   */
  double softLog(tcPDPtr p) const;

  /**
   * The Kbar^{qq} contribution
   */
  double KBarqq() const;

  /**
   * The Ktilde^{qq} contribution
   */
  double KTildeqq() const;

  /**
   * The P^{qq} contribution
   */
  double Pqq() const;

  /**
   * The Kbar^{qg} contribution
   */
  double KBarqg() const;

  /**
   * The Ktilde^{qg} contribution
   */
  double KTildeqg() const;

  /**
   * The P^{qg} contribution
   */
  double Pqg() const;

  /**
   * The Kbar^{gq} contribution
   */
  double KBargq() const;

  /**
   * The Ktilde^{gq} contribution
   */
  double KTildegq() const;

  /**
   * The P^{gq} contribution
   */
  double Pgq() const;

  /**
   * The Kbar^{gg} contribution
   */
  double KBargg() const;

  /**
   * The Ktilde^{gg} contribution
   */
  double KTildegg() const;

  /**
   * The P^{gg} contribution
   */
  double Pgg() const;

  /**
   * Get all contributions for the indicated incoming parton.
   */
  double sumParton(int id) const;

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
   * K_q
   */
  mutable double KQuark;

  /**
   * K_g
   */
  mutable  double KGluon;

  /**
   * The scale to be used.
   */
  mutable Energy2 scale;

  /**
   * The PDF to be used.
   */
  mutable tcPDFPtr pdf;

  /**
   * The incoming parton type considered.
   */
  mutable tcPDPtr particle;

  /**
   * The x to be used.
   */
  mutable double x;

  /**
   * The convolution variable to be used.
   */
  mutable double z;

  /**
   * Cache PDFs evaluated at x and x/z
   */
  mutable map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> > pdfCache;

  /**
   * The currently considered incoming parton.
   */
  mutable tcPDPtr parton;

  /**
   * Get a PDF value at x
   */
  double PDFx(tcPDPtr) const;

  /**
   * Get a PDF value at x/z
   */
  double PDFxByz(tcPDPtr) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipolePKOperator> initDipolePKOperator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipolePKOperator & operator=(const DipolePKOperator &);

};

}

#endif /* HERWIG_DipolePKOperator_H */
