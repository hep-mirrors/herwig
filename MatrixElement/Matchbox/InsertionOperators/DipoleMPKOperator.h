// -*- C++ -*-
//
// DipoleMPKOperator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleMPKOperator_H
#define HERWIG_DipoleMPKOperator_H
//
// This is the declaration of the DipoleMPKOperator class.
//

#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "ThePEG/PDF/PDF.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Daniel Rauch, Christian Reuschle,
 *         Johannes Bellm
 *
 * \brief DipoleMPKOperator implements the P+K
 * insertion operator for the massive case.
 * DipoleMPKOperator does not apply for dimen-
 * sional reduction.
 *
 */
class DipoleMPKOperator: public MatchboxInsertionOperator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleMPKOperator();

  /**
   * The destructor.
   */
  virtual ~DipoleMPKOperator();
  //@}

public:

  /**
   * Set the XComb object steering the Born matrix
   * element this class represents virtual corrections to.
   */
  virtual void setXComb(tStdXCombPtr xc);

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
   * Return true, if contributions exist to
   * the given parton.
   */
  bool applyNotMassless(tcPDPtr) const;

  /**
   * Return true, if this virtual correction
   * applies to the given process.
   */
  virtual bool apply(const cPDVector&) const;

  /**
   * The initial-final contribution
   * [1/(1-z)]_+ + \delta(1-z)
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

  //////////////////////////////////////

  /**
   * Kscript^{a,a'}_j terms with j=quark
   */

  /**
   * [J^a_{gQ}(z,\mu_Q^2)]_+
   * a,a' = quark for later
   * use of this function
   * in Kscript^{qq}_q
   */
  double Ja_gQplus(double muQ2) const;

  /**
   * [1/(1-z)]_+ * log( (2-z)/(2-z+\mu_Q^2) )
   * a,a' = quark for later use of this function
   * in Kscript^{qq}_q
   */
  double gammaSoft2(double muQ2) const;

  /**
   * The Kscript^{qq}_q contribution
   */
  double Kscriptqq_q(Energy2 sja, Energy2 mj2) const;

  /**
   * The Kscript^{qg}_q contribution
   */
  double Kscriptqg_q(Energy2 sja, Energy2 mj2) const;

  /**
   * The Kscript^{gq}_q contribution
   */
  double Kscriptgq_q() const;

  /**
   * The Kscript^{gg}_q contribution
   */
  double Kscriptgg_q(Energy2 sja, Energy2 mj2) const;

  /**
   * The Kscriptbar^{qq}_q contribution (B.17)
   */
  double Kscriptbarqq_q(Energy2 Qja2, Energy2 mj2) const;

  /**
   * The Kscriptbar^{qg}_q contribution (B.17)
   */
  double Kscriptbarqg_q(Energy2 Qja2, Energy2 mj2) const;

  /**
   * The Kscriptbar^{gq}_q contribution (B.17)
   */
  double Kscriptbargq_q() const;

  /**
   * The Kscriptbar^{gg}_q contribution (B.17)
   */
  double Kscriptbargg_q(Energy2 Qja2, Energy2 mj2) const;

  ////////////////////////////

  /**
   * Kscript^{a,a'}_j terms with j=gluon
   * The ones for a!=a' will return zero
   */

  /**
   * J^{a;NS}_{Q\bar{Q}}(\mu_Q^2)
   * Not folded with 1/z*PDF(x/z)*\Theta(z-x)
   */
  double JaNS_QQ(double muQ2) const;

  /**
   * [J^a_{Q\bar{Q}}(z,\mu_Q^2)]_{z_+}
   */
  double Ja_QQzplus(double muQ2, int F, double zplus) const;

  /**
   * The Kscript^{qq}_g contribution
   */
  // double Kscriptqq_g(Energy2 sja) const;
  double Kscriptqq_g(Energy2 sja, double lambda) const;

  /**
   * The Kscript^{qg}_g contribution
   */
  double Kscriptqg_g() const;

  /**
   * The Kscript^{gq}_g contribution
   */
  double Kscriptgq_g() const;

  /**
   * The Kscript^{gg}_g contribution
   * equals the Kscript^{qq}_g contribution
   */
  // double Kscriptgg_g(Energy2 sja) const;
  double Kscriptgg_g(Energy2 sja, double lambda) const;

  /**
   * The Kscriptbar^{qq}_g contribution (B.18)
   */
  // double Kscriptbarqq_g(Energy2 Qja2) const;
  double Kscriptbarqq_g(Energy2 Qja2, double lambda) const;

  /**
   * The Kscriptbar^{qg}_g contribution (B.18)
   */
  double Kscriptbarqg_g() const;

  /**
   * The Kscriptbar^{gq}_g contribution (B.18)
   */
  double Kscriptbargq_g() const;

  /**
   * The Kscriptbar^{gg}_g contribution (B.18)
   */
  // double Kscriptbargg_g(Energy2 Qja2) const;
  double Kscriptbargg_g(Energy2 Qja2, double lambda) const;

  //////////////////////////////////////

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
  double KQuark;

  /**
   * K_g
   */
  double KGluon;

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
   * Cache PDFs evaluated at x, x/z and x/z_+,
   * z_+ depends hereby on the masses of the heavy flavours,
   * so, in contrast to the massless case, we need a vector
   * of size 2+n_F(heavy) instead of the pair<double,double>
   */
  mutable map<pair<tcPDFPtr,tcPDPtr>,vector<double> > pdfCache;

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

  //////////////////////////////////////

  /**
   * Get a PDF value at x/z_+
   */
  double PDFxByzplus(tcPDPtr,int,double) const;

  //////////////////////////////////////

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleMPKOperator> initDipoleMPKOperator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleMPKOperator & operator=(const DipoleMPKOperator &);

};

}

#endif /* HERWIG_DipoleMPKOperator_H */
