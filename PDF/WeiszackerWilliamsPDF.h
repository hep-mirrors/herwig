// -*- C++ -*-
#ifndef HERWIG_WeiszackerWilliamsPDF_H
#define HERWIG_WeiszackerWilliamsPDF_H
//
// This is the declaration of the WeiszackerWilliamsPDF class.
//

#include "ThePEG/PDF/PDFBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the WeiszackerWilliamsPDF class.
 *
 * @see \ref WeiszackerWilliamsPDFInterfaces "The interfaces"
 * defined for WeiszackerWilliamsPDF.
 */
class WeiszackerWilliamsPDF: public PDFBase {

public:

  /**
   *  Default constructor
   */
  WeiszackerWilliamsPDF();

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
   * logarithmic momentum fraction \a l \f$(l=\log(1/x)\$f. The \a
   * particle is assumed to have a virtuality \a particleScale.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\$f. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * will return zero.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * Generate scale (as a fraction of the maximum scale). If the PDF
   * contains strange peaks which can be difficult to handle, this
   * function may be overwritten to return an appropriate scale
   * \f$Q^2/Q^2_{\max}\f$ for a \a z uniformly distributed in
   * ]0,1[. Also the jacobobian of the \f$Q^2/Q^2_{\max}\rightarrow
   * z\f$ variable transformation must multiply the \a jacobian
   * argument. The default version will simply use the function
   * \f$Q^2/Q^2_{\max} = (Q^2_{\max}/Q^2_{\min})^(z-1)\f$ or, if
   * \f$Q^2_{\min}\f$ is zero, \f$Q^2/Q^2_{\max} = z\f$ (where the
   * limits are set by \a cut).
   */
  virtual double flattenScale(tcPDPtr particle, tcPDPtr parton,
			       const PDFCuts & cut, double l, double z,
			       double & jacobian) const;

  /**
   * Generate a momentum fraction. If the PDF contains strange peaks
   * which can be difficult to handle, this function may be
   * overwritten to return an appropriate \f$l=\log(1/x)\f$ for a \a z
   * uniformly distributed in ]0,1[. Also the jacobobian of the
   * \f$l\rightarrow z\f$ variable transformation must in the function
   * multiply the \a jacobian argument. The default version will
   * simply use the function \f$l(z) = l_{\min} +
   * z*(l_{\max}-l_{\min})\f$ (where the limits are set by \a cut).
   */
  virtual double flattenL(tcPDPtr particle, tcPDPtr parton, const PDFCuts &cut,
			  double z, double & jacobian) const;
  //@}


public:

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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<WeiszackerWilliamsPDF> initWeiszackerWilliamsPDF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeiszackerWilliamsPDF & operator=(const WeiszackerWilliamsPDF &);

private:

  /**
   *  Minimum \f$Q^2\f$ for the photon
   */
  Energy2 _q2min;

  /**
   *  Maximum \f$Q^2\f$ for the photon
   */
  Energy2 _q2max;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of WeiszackerWilliamsPDF. */
template <>
struct BaseClassTrait<Herwig::WeiszackerWilliamsPDF,1> {
  /** Typedef of the first base class of WeiszackerWilliamsPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the WeiszackerWilliamsPDF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::WeiszackerWilliamsPDF>
  : public ClassTraitsBase<Herwig::WeiszackerWilliamsPDF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::WeiszackerWilliamsPDF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * WeiszackerWilliamsPDF is implemented. It may also include several, space-separated,
   * libraries if the class WeiszackerWilliamsPDF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLeptonPDF.so"; }
};

/** @endcond */

}

#endif /* HERWIG_WeiszackerWilliamsPDF_H */
