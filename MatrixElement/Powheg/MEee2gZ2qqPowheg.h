// -*- C++ -*-
#ifndef HERWIG_MEee2gZ2qqPowheg_H
#define HERWIG_MEee2gZ2qqPowheg_H
//
// This is the declaration of the MEee2gZ2qqPowheg class.
//

#include "Herwig++/MatrixElement/Lepton/MEee2gZ2qq.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEee2gZ2qqPowheg class.
 *
 * @see \ref MEee2gZ2qqPowhegInterfaces "The interfaces"
 * defined for MEee2gZ2qqPowheg.
 */
class MEee2gZ2qqPowheg: public MEee2gZ2qq {

public:

  /**
   * The default constructor.
   */
  MEee2gZ2qqPowheg() : contrib_(1), zPow_(0.5), yPow_(0.9)
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;
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
  static ClassDescription<MEee2gZ2qqPowheg> initMEee2gZ2qqPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2gZ2qqPowheg & operator=(const MEee2gZ2qqPowheg &);

private:

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Phase-space sampling for z
   */
  double zPow_;

  /**
   *  Phase-space sampling for y
   */
  double yPow_;

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double z_;

  /**
   *  The \f$y\f$ angular variable
   */
  double y_;

  /**
   *  The azimuth
   */
  double phi_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2gZ2qqPowheg. */
template <>
struct BaseClassTrait<Herwig::MEee2gZ2qqPowheg,1> {
  /** Typedef of the first base class of MEee2gZ2qqPowheg. */
  typedef Herwig::MEee2gZ2qq NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2gZ2qqPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2gZ2qqPowheg>
  : public ClassTraitsBase<Herwig::MEee2gZ2qqPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2gZ2qqPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2gZ2qqPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEee2gZ2qqPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so HwPowhegMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2gZ2qqPowheg_H */
