// -*- C++ -*-
#ifndef HERWIG_MPIPDF_H
#define HERWIG_MPIPDF_H
//
// This is the declaration of the MPIPDF class.
//

#include "ThePEG/PDF/PDFBase.h"
#include "MPIPDF.fh"

namespace Herwig {
using namespace ThePEG;
/**
 * Here is the documentation of the MPIPDF class. It defines
 * a modified pdf which uses an existing pdf object to add
 * modifications like removing the valence part of it, which
 * is needed in the backward evolution of secondary scatters.
 *
 * @see \ref MPIPDFInterfaces "The interfaces"
 * defined for MPIPDF.
 */
class MPIPDF: public PDFBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor, which shouldn't be called
   */
  inline MPIPDF();

  /**
   * The constructor which takes a PDF object as argument, to work with.
   */
  inline MPIPDF(tcPDFPtr orig);

  /**
   * The copy constructor.
   */
  inline MPIPDF(const MPIPDF &);

  /**
   * The destructor.
   */
  virtual ~MPIPDF();
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MPIPDF> initMPIPDF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MPIPDF & operator=(const MPIPDF &);

  /**
   * pointer to the underlying ThePEG::PDFBase object, we are modifying.
   */
  tcPDFPtr thePDF;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MPIPDF. */
template <>
struct BaseClassTrait<Herwig::MPIPDF,1> {
  /** Typedef of the first base class of MPIPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MPIPDF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MPIPDF>
  : public ClassTraitsBase<Herwig::MPIPDF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MPIPDF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MPIPDF is implemented. It may also include several, space-separated,
   * libraries if the class MPIPDF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPIPDF.so"; }
};

/** @endcond */

}

#include "MPIPDF.icc"
#ifndef HERWIG_TEMPLATES_IN_CC_FILE
// #include "MPIPDF.tcc"
#endif

#endif /* HERWIG_MPIPDF_H */
