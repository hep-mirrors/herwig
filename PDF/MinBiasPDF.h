// -*- C++ -*-
//
// MinBiasPDF.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MinBiasPDF_H
#define HERWIG_MinBiasPDF_H
//
// This is the declaration of the MinBiasPDF class.
//

#include "ThePEG/PDF/PDFBase.h"
#include "MinBiasPDF.fh"

namespace Herwig {
using namespace ThePEG;
/**
 * The MinBiasPDF class defines
 * a modified pdf which uses an existing pdf object to add
 * modifications like removing the valence part of it, which
 * is needed in the backward evolution of secondary scatters.
 *
 * \author Manuel B\"ahr
 *
 * @see \ref MinBiasPDFInterfaces "The interfaces"
 * defined for MinBiasPDF.
 */
class MinBiasPDF: public PDFBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MinBiasPDF() : thePDF(PDFPtr()) {}

  /**
   * The copy constructor.
   */
  MinBiasPDF(const MinBiasPDF & x) : 
    PDFBase(x), thePDF(x.thePDF) {}

  /**
   * The destructor.
   */
  virtual ~MinBiasPDF() {}

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

  /** return the underlying PDFBase pointer*/
  tcPDFPtr originalPDF() const {return thePDF;}

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
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MinBiasPDF> initMinBiasPDF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MinBiasPDF & operator=(const MinBiasPDF &);

  /**
   * pointer to the underlying ThePEG::PDFBase object, we are modifying.
   */
  PDFPtr thePDF;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MinBiasPDF. */
template <>
struct BaseClassTrait<Herwig::MinBiasPDF,1> {
  /** Typedef of the first base class of MinBiasPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MinBiasPDF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MinBiasPDF>
  : public ClassTraitsBase<Herwig::MinBiasPDF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MinBiasPDF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MinBiasPDF is implemented. It may also include several, space-separated,
   * libraries if the class MinBiasPDF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MinBiasPDF_H */
