// -*- C++ -*-
//
// PDFRatio.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PDFRatio_H
#define HERWIG_PDFRatio_H
//
// This is the declaration of the PDFRatio class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/PDF/PDF.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 * 
 * \brief PDFRatio implements numerically stable PDF ratios.
 *
 * @see \ref PDFRatioInterfaces "The interfaces"
 * defined for PDFRatio.
 */
class PDFRatio: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PDFRatio();

  /**
   * The destructor.
   */
  virtual ~PDFRatio();
  //@}

public:

  /**
   * For the given PDF, scale and partons from and to and
   * x,z values return the ratio xf_to(x/z) / xf_from(x)
   */
  double operator() (const PDF& pdf,
		     Energy2 scale,
		     tcPDPtr from, tcPDPtr to,
		     double x, double z) const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The x from which on extrapolation should
   * be done for valence partons.
   */
  double theValenceExtrapolation;

  /**
   * The x from which on extrapolation should
   * be done for sea partons.
   */
  double theSeaExtrapolation;

  /**
   * The scale below which the PDF will be frozen
   */
  Energy theFreezingScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PDFRatio> initPDFRatio;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PDFRatio & operator=(const PDFRatio &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PDFRatio. */
template <>
struct BaseClassTrait<Herwig::PDFRatio,1> {
  /** Typedef of the first base class of PDFRatio. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PDFRatio class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PDFRatio>
  : public ClassTraitsBase<Herwig::PDFRatio> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PDFRatio"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PDFRatio is implemented. It may also include several, space-separated,
   * libraries if the class PDFRatio depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PDFRatio_H */
