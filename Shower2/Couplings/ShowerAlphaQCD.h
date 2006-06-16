// -*- C++ -*-
#ifndef THEPEG_ShowerAlphaQCD_H
#define THEPEG_ShowerAlphaQCD_H
//
// This is the declaration of the ShowerAlphaQCD class.
//

#include "ShowerAlpha.h"
#include "ThePEG/Config/Constants.h"
#include "ShowerAlphaQCD.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value() and overestimateValue() for the 
 *  strong coupling.
 *
 *  A  number of different options for the running of the coupling
 *  and its initial definition are supported.
 *
 * @see \ref ShowerAlphaQCDInterfaces "The interfaces"
 * defined for ShowerAlphaQCD.
 */
class ShowerAlphaQCD: public ShowerAlpha {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerAlphaQCD();
  //@}

public:

  /**
   *  Methods to return the coupling
   */
  //@{
  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   * @param scale The scale
   * @return The coupling
   */
  virtual double value(const Energy2 scale) const;

  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double overestimateValue() const;

  /**
   *  Return the ratio of the coupling at the scale to the overestimated value
   */
  virtual double ratio(const Energy2 scale) const;
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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   *  Member functions which calculate the coupling
   */
  //@{
  /**
   * The 1,2,3-loop parametrization of \f$\alpha_S\f$.
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
  inline double alphaS(Energy q, Energy lam, short nf) const; 

  /**
   * The derivative of \f$\alpha_S\f$ with respect to \f$\ln(Q^2/\Lambda^2)\f$
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
  inline double derivativealphaS(Energy q, Energy lam, short nf) const; 

  /**
   * Compute the value of \f$Lambda\f$ needed to get the input value of
   * the strong coupling at the scale given for the given number of flavours
   * using the Newton-Raphson method
   * @param match The scale for the coupling
   * @param alpha The input coupling
   * @param nflav The number of flavours
   */
  inline Energy computeLambda(Energy match,double alpha,unsigned int nflav) const;

  /**
   * Return the value of \f$\Lambda\f$ and the number of flavours at the scale.
   * @param q The scale
   * @return The number of flavours at the scale and \f$\Lambda\f$.
   */
  inline pair<short, Energy> getLamNfTwoLoop(Energy q) const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerAlphaQCD> initShowerAlphaQCD;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlphaQCD & operator=(const ShowerAlphaQCD &);

private:

  /**
   *  Minimum value of the scale
   */
  Energy _qmin;

  /**
   *  Parameter controlling the behaviour of \f$\alpha_S\f$ in the non-perturbative
   *  region.
   */ 
  int _asType;

  /**
   *  Thresholds for the different number of flavours 
   */
  vector<Energy> _thresholds;

  /**
   *  \f$\Lambda\f$ for the different number of flavours
   */
  vector<Energy> _lambda;

  /**
   *  Option for the number of loops
   */
  unsigned int _nloop;

  /**
   *  Option for the translation between \f$\Lambda_{\bar{MS}}\f$ and
   *  \f$\Lambda_{\rm Herwig}\f$
   */
  bool _lambdaopt;

  /**
   *  Input value of Lambda
   */
  Energy _lambdain;

  /**
   *  Input value of \f$alpha_S(M_Z)\f$
   */
  double _alphain;

  /**
   *  Option for the calculation of Lambda from input parameters
   */
  bool _inopt;

  /**
   *  Tolerance for discontinuities at the thresholds
   */
  double _tolerance;

  /**
   *  Maximum number of iterations for the Newton-Raphson method to converge
   */
  unsigned int _maxtry;

  /**
   *  The minimum value of the coupling
   */
  double _alphamin;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerAlphaQCD. */
template <>
struct BaseClassTrait<Herwig::ShowerAlphaQCD,1> {
  /** Typedef of the first base class of ShowerAlphaQCD. */
  typedef Herwig::ShowerAlpha NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerAlphaQCD class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerAlphaQCD>
  : public ClassTraitsBase<Herwig::ShowerAlphaQCD> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ShowerAlphaQCD"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerAlphaQCD is implemented. It may also include several, space-separated,
   * libraries if the class ShowerAlphaQCD depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "ShowerAlphaQCD.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQCD.tcc"
#endif

#endif /* THEPEG_ShowerAlphaQCD_H */
