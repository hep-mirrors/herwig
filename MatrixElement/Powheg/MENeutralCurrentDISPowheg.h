// -*- C++ -*-
#ifndef HERWIG_MENeutralCurrentDISPowheg_H
#define HERWIG_MENeutralCurrentDISPowheg_H
//
// This is the declaration of the MENeutralCurrentDISPowheg class.
//

#include "Herwig++/MatrixElement/DIS/MENeutralCurrentDIS.h"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Typedef for BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;
/**
 * Here is the documentation of the MENeutralCurrentDISPowheg class.
 *
 * @see \ref MENeutralCurrentDISPowhegInterfaces "The interfaces"
 * defined for MENeutralCurrentDISPowheg.
 */
class MENeutralCurrentDISPowheg: public MENeutralCurrentDIS {

public:

  /**
   * The default constructor.
   */
  MENeutralCurrentDISPowheg();

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

  /**
   *  The NLO weight
   */
  double NLOWeight() const;

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
  static ClassDescription<MENeutralCurrentDISPowheg> initMENeutralCurrentDISPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MENeutralCurrentDISPowheg & operator=(const MENeutralCurrentDISPowheg &);

private:

  /**
   *  The Born variables
   */
  //@{
  /**
   *  \f$x_B\f$
   */
  double _xB;

  /**
   *   \f$Q^2\f$
   */
  Energy2 _q2;
  //@}

  /**
   *  The radiative variables
   */
  //@{
  /**
   *  The \f$x_p\f$ or \f$z\f$ real integration variable
   */
  double _xp;
  //@}

  /**
   *  The hadron
   */
  tcBeamPtr _hadron;

  /**
   * Selects a dynamic or fixed factorization scale
   */
  unsigned int scaleOpt_;

  /**
   * The factorization scale 
   */
  Energy muF_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Jacobian for \f$x_p\f$ integral
   */
  double jac_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MENeutralCurrentDISPowheg. */
template <>
struct BaseClassTrait<Herwig::MENeutralCurrentDISPowheg,1> {
  /** Typedef of the first base class of MENeutralCurrentDISPowheg. */
  typedef Herwig::MENeutralCurrentDIS NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MENeutralCurrentDISPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MENeutralCurrentDISPowheg>
  : public ClassTraitsBase<Herwig::MENeutralCurrentDISPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MENeutralCurrentDISPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MENeutralCurrentDISPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MENeutralCurrentDISPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEDIS.so HwPowhegMEDIS.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MENeutralCurrentDISPowheg_H */
