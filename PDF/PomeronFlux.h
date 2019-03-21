// -*- C++ -*-
#ifndef Herwig_PomeronFlux_H
#define Herwig_PomeronFlux_H
//
// This is the declaration of the PomeronFlux class.
//

#include "ThePEG/PDF/PDFBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the PomeronFlux class.
 *
 * @see \ref PomeronFluxInterfaces "The interfaces"
 * defined for PomeronFlux.
 */
class PomeronFlux: public PDFBase {

public:

  /**
   *  Default constructor
   */
  PomeronFlux();

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
		     double l, Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\f$. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * will return zero.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = ZERO) const;

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

  /**
   * Set normalization constant Ap of the pomeron flux. 
   */
   virtual void doinit();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static ClassDescription<PomeronFlux> initPomeronFlux;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PomeronFlux & operator=(const PomeronFlux &) = delete;

private:
  
  /**
   *  Integrated pomeron/regeon flux over qq in the range qqmin ... qqmax 
   *  (with normalization constant Ap = 1 GeV2). 
   *  The flux normalization constant is obtained \f$x_p*intFx_p = 1\f$, 
   *  where x_p = 0.003.
   */
  Energy2 intxFx(double x, Energy2 qqmin, Energy2 qqmax, 
		 double alfa0, InvEnergy2 alfap,  InvEnergy2 beta) const;

  /**
   * Helper function in doinit() which sets the parameters of the 
   * pomeron/regeon flux according to user setup.
   */
  void setFluxPar();
  

private:

  /**
   *  Minimum \f$Q^2\f$ for the pomeron/reggeon
   */
  Energy2 q2min_;

  /**
   *  Maximum \f$Q^2\f$ for the pomeron/reggeon
   */
  Energy2 q2max_;

  /**
   * Cut on the minimum xi
   */
  double xiMin_;
  
  /**
   * Cut on the maximum xi
   */
  double xiMax_;
  
  /**
   * Pomeron intercept  
   */
  double alfa0P_; 

  /**
   *  Pomeron slope 
   */
  InvEnergy2 alfapP_; 

  /**
   * Parameter of pomeron flux Bp    
   */
  InvEnergy2 betaP_;
  
  /**
   * Normalization constant of the pomeron flux    
   */
  InvEnergy2 normP_;

  /**
   * Reggeon intercept  
   */
  double alfa0R_; 

  /**
   *  Reggeon slope 
   */
  InvEnergy2 alfapR_; 

  /**
   * Parameter of reggeon flux Bp    
   */
  InvEnergy2 betaR_;
  
  /**
   * Normalization constant of the reggeon flux    
   */
  InvEnergy2 normR_;
  
  /**
   * Factor of the pomeron flux    
   */
  double nR_;
  
  /**
   * Switch between pomeron/regeon structure function fits.
   * The flux parameters are set according to the fit choice.     
   */
  int PDFFit_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PomeronFlux. */
template <>
struct BaseClassTrait<Herwig::PomeronFlux,1> {
  /** Typedef of the first base class of PomeronFlux. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PomeronFlux class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PomeronFlux>
  : public ClassTraitsBase<Herwig::PomeronFlux> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PomeronFlux"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PomeronFlux is implemented. It may also include several, space-separated,
   * libraries if the class PomeronFlux depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPomeronFlux.so"; }
};

/** @endcond */

}

#endif /* Herwig_PomeronFlux_H */
