// -*- C++ -*-
#ifndef HERWIG_MEPP2HiggsVBFPowheg_H
#define HERWIG_MEPP2HiggsVBFPowheg_H
//
// This is the declaration of the MEPP2HiggsVBFPowheg class.
//

#include "Herwig/MatrixElement/Hadron/MEPP2HiggsVBF.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2HiggsVBFPowheg class.
 *
 * @see \ref MEPP2HiggsVBFPowhegInterfaces "The interfaces"
 * defined for MEPP2HiggsVBFPowheg.
 */
class MEPP2HiggsVBFPowheg: public Herwig::MEPP2HiggsVBF {

public:

  /**
   * The default constructor.
   */
  MEPP2HiggsVBFPowheg();

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

  /**
   *  Leading order matrix element
   */
  Energy4 loMatrixElement(const Lorentz5Momentum &p1,
			  const Lorentz5Momentum &p2,
			  const Lorentz5Momentum &q1,
			  const Lorentz5Momentum &q2,
			  double G1, double G2) const;

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
  static ClassDescription<MEPP2HiggsVBFPowheg> initMEPP2HiggsVBFPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsVBFPowheg & operator=(const MEPP2HiggsVBFPowheg &);

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
   *  Partons
   */
  tcPDPtr _partons[5];

  mutable Energy2 _q2;
  //@}

  /**
   *  The radiative variables
   */
  //@{
  /**
   *  The \f$x_p\f$ or \f$z\f$ real integration variable
   */
  double  _xp;

  /**
   *  The \f$z_p\f$ real integration variable 
   */
  double _zp;

  /**
   *  The \f$fi\f$ real integration variable 
   */
  double _phi;
  //@}

  /**
   *  The variables to get the right boost
   */
  //@{
  /**
   *  LO momenta 
   */
  Lorentz5Momentum _loMomenta[4];
  /**
   *  The transfered (virtual boson) momentum 
   */
  mutable Lorentz5Momentum _pa;

  /**
   *  The incoming quark momentum 
   */
  mutable Lorentz5Momentum _pb;

  /**
   *  The outgoing quark momentum
   */
  Lorentz5Momentum _pc;

  /**
   *  The incoming quark momentum 
   */
  Lorentz5Momentum _pbother;

  /**
   *  The outgoing quark momentum
   */
  Lorentz5Momentum _pcother;
  //@}

  /**
   *  Electroweak parameters
   */
  //@{
  /**
   *  The square of the Z mass
   */
  Energy2 _mz2;

  /**
   *  The square of the W mass
   */
  Energy2 _mw2;
  //@}

  /**
   *  The first hadron
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
 *  base classes of MEPP2HiggsVBFPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2HiggsVBFPowheg,1> {
  /** Typedef of the first base class of MEPP2HiggsVBFPowheg. */
  typedef Herwig::MEPP2HiggsVBF NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2HiggsVBFPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2HiggsVBFPowheg>
  : public ClassTraitsBase<Herwig::MEPP2HiggsVBFPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2HiggsVBFPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2HiggsVBFPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2HiggsVBFPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2HiggsVBFPowheg_H */
