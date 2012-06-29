// -*- C++ -*-
#ifndef HERWIG_FortranAlphaQCD_H
#define HERWIG_FortranAlphaQCD_H
//
// This is the declaration of the FortranAlphaQCD class.
//

#include "ShowerAlpha.h"
#include "ThePEG/Config/Constants.h"
#include "FortranAlphaQCD.fh"

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
 * @see \ref FortranAlphaQCDInterfaces "The interfaces"
 * defined for FortranAlphaQCD.
 */
class FortranAlphaQCD: public ShowerAlpha {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline FortranAlphaQCD();
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

  /**
   *  Get the value of \f$\Lambda_{\rm QCD}\f$ for three flavours
   */
  inline Energy lambdaQCDThree() const;

  /**
   *  Get the value of \f$\Lambda_{\rm QCD}\f$ for five flavours
   */
  inline Energy lambdaQCDFive() const;

  /**
   *  The one-loop coupling with five flavour \f$\beta\f$ and three flavour
   *  \f$\Lambda_{\rm QCD}\f$
   */
  double oneLoopValue(const Energy2 scale) const;

  /**
   *  The ratio of the two-loop to one-loop results
   */
  double twoLoopRatio(const Energy2 scale) const;

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
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FortranAlphaQCD> initFortranAlphaQCD;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranAlphaQCD & operator=(const FortranAlphaQCD &);

private:

  /**
   *  Thresholds for the different number of flavours 
   */
  vector<Energy> _thresholds;

  /**
   *  Input value of \f$\Lambda\f$
   */
  Energy _lambdain;

  /**
   *  5-flavour \f$\Lambda\f$
   */
  Energy _lambda5;

  /**
   *  3-flavour \f$\Lambda\f$
   */
  Energy _lambda3;

  /**
   *  Maximum value
   */
  double _alphamax;

  /**
   *  Maximum number of iterations for the Newton-Raphson method to converge
   */
  unsigned int _maxtry;

  /**
   *  Matching coefficients
   */
  vector<double> _match;

  /**
   *  First \f$\beta\f$ function coefficient
   */
  vector<double> _bcoeff;

  /**
   *  Second \f$\beta\f$ function coefficient
   */
  vector<double> _ccoeff;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FortranAlphaQCD. */
template <>
struct BaseClassTrait<Herwig::FortranAlphaQCD,1> {
  /** Typedef of the first base class of FortranAlphaQCD. */
  typedef Herwig::ShowerAlpha NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FortranAlphaQCD class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FortranAlphaQCD>
  : public ClassTraitsBase<Herwig::FortranAlphaQCD> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FortranAlphaQCD"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FortranAlphaQCD is implemented. It may also include several, space-separated,
   * libraries if the class FortranAlphaQCD depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "FortranAlphaQCD.icc"

#endif /* HERWIG_FortranAlphaQCD_H */
