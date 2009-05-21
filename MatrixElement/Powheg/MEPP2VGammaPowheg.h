// -*- C++ -*-
#ifndef HERWIG_MEPP2VGammaPowheg_H
#define HERWIG_MEPP2VGammaPowheg_H
//
// This is the declaration of the MEPP2VGammaPowheg class.
//

#include "Herwig++/MatrixElement/Hadron/MEPP2VGamma.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2VGammaPowheg class implements the next-to-leading
 * order matrix elements for $q\bar q \to W^\pm/Z^0\gamma\f$
 * in the Powheg scheme.
 *
 * @see \ref MEPP2VGammaPowhegInterfaces "The interfaces"
 * defined for MEPP2VGammaPowheg.
 */
class MEPP2VGammaPowheg: public MEPP2VGamma {

public:

  /**
   * The default constructor.
   */
  MEPP2VGammaPowheg();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2VGammaPowheg> initMEPP2VGammaPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VGammaPowheg & operator=(const MEPP2VGammaPowheg &);

protected:

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

private:

  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int _contrib;
  //@}

  /**
   *  Choice of the scale
   */
  //@{
  /**
   *  Type of scale
   */
  unsigned int _scaleopt;

  /**
   *  Fixed scale if used
   */
  Energy _fixedScale;

  /**
   *  Prefactor if variable scale used
   */
  double _scaleFact;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VGammaPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2VGammaPowheg,1> {
  /** Typedef of the first base class of MEPP2VGammaPowheg. */
  typedef Herwig::MEPP2VGamma NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VGammaPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VGammaPowheg>
  : public ClassTraitsBase<Herwig::MEPP2VGammaPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VGammaPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VGammaPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VGammaPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VGammaPowheg_H */
