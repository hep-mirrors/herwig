// -*- C++ -*-
//
// SatPDF.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SatPDF_H
#define HERWIG_SatPDF_H
//
// This is the declaration of the SatPDF class.
//

#include "ThePEG/PDF/PDFBase.h"

namespace Herwig {
using namespace ThePEG;
/**
 * The SatPDF class defines
 * a modified pdf which uses an existing pdf object to add
 * modifications like removing the valence part of it, which
 * is needed in the backward evolution of secondary scatters.
 *
 * \author Manuel B\"ahr
 *
 * @see \ref SatPDFInterfaces "The interfaces"
 * defined for SatPDF.
 */
class SatPDF: public PDFBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SatPDF() : thePDF(PDFPtr()), theX0(1E-4), theExp(0.0) {}

  /**
   * The copy constructor.
   */
  inline SatPDF(const SatPDF & x) : 
    PDFBase(x), thePDF(x.thePDF), theX0(x.theX0), theExp(x.theExp) {}

  /**
   * The destructor.
   */
  virtual ~SatPDF();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of partons from
   * the given \a particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the partons which this PDF may extract from the given
   * \a particle.
   */
  virtual cPDVector partons(tcPDPtr particle) const;

  /**
   * The density. Return the pdf for the given \a parton inside the
   * given \a particle for the virtuality \a partonScale and momentum
   * fraction \a x. The \a particle is assumed to have a virtuality \a
   * particleScale.
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps=0.0, Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and momentum fraction \a x. The \a particle is
   * assumed to have a virtuality \a particleScale. If not overidden
   * by a sub class this implementation will assume that the
   * difference between a quark and anti-quark distribution is due do
   * valense quarks, but return zero for anything else.
   */
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps=0.0, Energy2 particleScale = ZERO) const;
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
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SatPDF & operator=(const SatPDF &) = delete;

  /**
   * pointer to the underlying ThePEG::PDFBase object, we are modifying.
   */
  PDFPtr thePDF;

  /**
   * x from where the extrapolation f(x) = f(theX0) * (x/theX0)**theExp
   * is used for the pdf's.
   */
  double theX0;

  /**
   * the exponent of the pdf extrapolation for small x (x < theX0).
   */
  double theExp;

};

}

#ifndef HERWIG_TEMPLATES_IN_CC_FILE
// #include "SatPDF.tcc"
#endif

#endif /* HERWIG_SatPDF_H */
