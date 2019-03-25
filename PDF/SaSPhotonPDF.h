// -*- C++ -*-
#ifndef Herwig_SaSPhotonPDF_H
#define Herwig_SaSPhotonPDF_H
//
// This is the declaration of the SaSPhotonPDF class.
//

#include "ThePEG/PDF/PDFBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SaSPhotonPDF class provides an interface to the 
 *
 * @see \ref SaSPhotonPDFInterfaces "The interfaces"
 * defined for SaSPhotonPDF.
 */
class SaSPhotonPDF: public PDFBase {

public:
  
  /**
   * The default constructor.
   */
  SaSPhotonPDF() : iset_(2), ip_(0) {}

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
   * given \a particle for the virtuality \a partonScale and
   * logarithmic momentum fraction \a l \f$(l=\log(1/x)\f$. The \a
   * particle is assumed to have a virtuality \a particleScale.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\f$. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * will return zero.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SaSPhotonPDF & operator=(const SaSPhotonPDF &) = delete;

private:

  /**
   *  PDF Set
   */
  mutable int iset_;

  /**
   *  Scheme used to evaluate off-shell anomalous component.
   */
  mutable int ip_;

};

}

#endif /* Herwig_SaSPhotonPDF_H */
