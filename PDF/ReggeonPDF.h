// -*- C++ -*-
//
// ReggeonPDF.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ReggeonPDF_H
#define HERWIG_ReggeonPDF_H

#include <ThePEG/PDF/PDFBase.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace Herwig {

using namespace ThePEG;

/** \ingroup PDF
 *
 *  Implementation of the ReggeonPDF PDFs
 *
 * @see \ref ReggeonPDFInterfaces "The interfaces"
 * defined for ReggeonPDF
 * This class is wrapper of Reggeon structure function. Which is most likely 
 * simulated by pion structure function.
 */
class ReggeonPDF : public PDFBase {

public:

  /**
   *  Default constructor
   */
  ReggeonPDF() : particleID_(111) {}

  /** @name Virtual functions from PDFBase */
  //@{
  /**
   * Return true if this PDF can handle the extraction of parton from the
   * given particle ie. if the particle is a proton or neutron.
   * @param particle The particle
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the parton types which are described by these parton
   * densities.
   * @param p The particle
   */
  virtual cPDVector partons(tcPDPtr p) const;

  /**
   * Return x times the pdf for the given parameters
   * @param particle The beam particle
   * @param parton The parton for which to return the PDF.
   * @param partonScale The scale at which to evaluate the PDF.
   * @param x The momentum fraction
   * @param eps ??? an unknown parameter from ThePEG.
   * @param particleScale The scale for the particle
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                     double x, double eps = 0.0,
                     Energy2 particleScale = ZERO) const;

  /**
   * Return x times the valence pdf for the given parameters
   * @param particle The beam particle
   * @param parton The parton for which to return the PDF.
   * @param partonScale The scale at which to evaluate the PDF.
   * @param x The momentum fraction
   * @param eps ??? an unknown parameter from ThePEG.
   * @param particleScale The scale for the particle
   */
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double eps = 0.0,
                      Energy2 particleScale = ZERO) const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ReggeonPDF> initReggeonPDF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ReggeonPDF & operator=(const ReggeonPDF &);

  /**
   * Pointer to the concrete PDF reggeon structure function. 
   */
  PDFPtr ptrPDF_;

  /**
   *  PDG code for the particle
   */
  long particleID_;

  /**
   *  Pointer to the particle
   */
  PDPtr particle_;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ReggeonPDF. */
template <>
struct BaseClassTrait<Herwig::ReggeonPDF,1> {
  /** Typedef of the first base class of ReggeonPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ReggeonPDF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ReggeonPDF>: public ClassTraitsBase<Herwig::ReggeonPDF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ReggeonPDF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ReggeonPDF is implemented. It may also include several, space-separated,
   * libraries if the class ReggeonPDF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwReggeonPDF.so"; }
};

/** @endcond */

}

#endif
