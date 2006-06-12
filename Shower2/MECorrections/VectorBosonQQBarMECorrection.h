// -*- C++ -*-
#ifndef HERWIG_VectorBosonQQBarMECorrection_H
#define HERWIG_VectorBosonQQBarMECorrection_H
//
// This is the declaration of the VectorBosonQQBarMECorrection class.
//

#include "MECorrectionBase.h"
#include "VectorBosonQQBarMECorrection.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the VectorBosonQQBarMECorrection class.
 *
 * @see \ref VectorBosonQQBarMECorrectionInterfaces "The interfaces"
 * defined for VectorBosonQQBarMECorrection.
 */
class VectorBosonQQBarMECorrection: public MECorrectionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline VectorBosonQQBarMECorrection();

  /**
   * The copy constructor.
   */
  inline VectorBosonQQBarMECorrection(const VectorBosonQQBarMECorrection &);

  /**
   * The destructor.
   */
  virtual ~VectorBosonQQBarMECorrection();
  //@}

public:

  /**
   *  Members to override those in the base class and implemented 
   *  the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   */
  virtual bool canHandle(ShowerTreePtr);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

  /**
   *  Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br);
  //@}

private:

  /**
   *  Apply the hard matrix element
   */
  vector<Lorentz5Momentum> applyHard(const ParticleVector &p);

  /**
   *  Get the weight for hard emission
   */
  double getHard(double &, double &);

  /**
   *  Set the \f$\rho\f$ parameter
   */
  inline void setRho(double);

  /**
   *  Set the \f$\tilde{\kappa}\f$ parameters symmetrically 
   */
  inline void setKtildeSymm();

  /**
   * Set second f$\tilde{\kappa}\f$, given the first.
   */
  inline void setKtilde2();

  /**
   *  Translate the variables from \f$x_q,x_{\bar{q}}\f$ to \f$\tilde{\kappa},z\f$
   */
  //@{
  /**
   *  Calculate \f$z\f$.
   */
  inline double getZfromX(double, double);

  /**
   *  Calculate \f$\tilde{\kappa}\f$.
   */
  inline double getKfromX(double, double);
  //@}

  /**
   * Calculate \f$x_{q},x_{\bar{q}}\f$ from \f$\tilde{\kappa},z\f$.
   * @param kt \f$\tilde{\kappa}\f$
   * @param z \f$z\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  inline void getXXbar(double kt, double z, double & x, double & xbar);

  /**
   *  Soft weight
   */
  //@{
  /**
   *  Soft quark weight calculated from \f$x_{q},x_{\bar{q}}\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  double qWeight(double x, double xbar); 

  /**
   *  Soft antiquark weight calculated from \f$x_{q},x_{\bar{q}}\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  double qbarWeight(double x, double xbar);

  /**
   * Soft quark weight calculated from \f$\tilde{q},z\f$
   * @param qtilde  \f$\tilde{q}\f$
   * @param z \f$z\f$
   */
  double qWeightX(Energy qtilde, double z);

  /**
   * Soft antiquark weight calculated from \f$\tilde{q},z\f$
   * @param qtilde  \f$\tilde{q}\f$
   * @param z \f$z\f$
   */
  double qbarWeightX(Energy qtilde, double z);
  //@}
  /**
   * ????
   */
  inline double u(double);

  /**
   *  Vector and axial vector parts of the matrix element
   */
  //@{
  /**
   *  Vector part of the matrix element
   */
  double MEV(double, double);

  /**
   *  Axial vector part of the matrix element
   */
  double MEA(double, double);

  /**
   * The matrix element, given \f$x_1\f$, \f$x_2\f$.
   * @param x1 \f$x_1\f$
   * @param x2 \f$x_2\f$
   */
  inline double PS(double x1, double x2);
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<VectorBosonQQBarMECorrection> initVectorBosonQQBarMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorBosonQQBarMECorrection & operator=(const VectorBosonQQBarMECorrection &);

private:

  /**
   * CM energy 
   */
  Energy d_Q;

  /**
   *  Quark mass
   */
  Energy d_m;

  /**
   * The rho parameter 
   */
  double d_rho;

  /**
   * The v parameter
   */
  double d_v;

  /**
   * The initial kappa-tilde values for radiation from the quark
   */
  double d_kt1;

  /**
   * The initial kappa-tilde values for radiation from the antiquark
   */
  double d_kt2;

  /**
   *  Cut-off parameter
   */
  static const double EPS;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VectorBosonQQBarMECorrection. */
template <>
struct BaseClassTrait<Herwig::VectorBosonQQBarMECorrection,1> {
  /** Typedef of the first base class of VectorBosonQQBarMECorrection. */
  typedef Herwig::MECorrectionBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VectorBosonQQBarMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VectorBosonQQBarMECorrection>
  : public ClassTraitsBase<Herwig::VectorBosonQQBarMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::VectorBosonQQBarMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VectorBosonQQBarMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class VectorBosonQQBarMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "VectorBosonQQBarMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorBosonQQBarMECorrection.tcc"
#endif

#endif /* HERWIG_VectorBosonQQBarMECorrection_H */
