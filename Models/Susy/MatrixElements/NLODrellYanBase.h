// -*- C++ -*-
#ifndef HERWIG_NLODrellYanBase_H
#define HERWIG_NLODrellYanBase_H
//
// This is the declaration of the NLODrellYanBase class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The NLODrellYanBase class provides a base class for the implementation
 * of Drell-Yan type processes in the POWHEG approach using Catani-Seymour
 * dipoles.
 *
 * @see \ref NLODrellYanBaseInterfaces "The interfaces"
 * defined for NLODrellYanBase.
 */
class NLODrellYanBase: public HwME2to2Base {

protected:

  /**
   *  Struct to return singular numbers
   */
  struct Singular {
    /**
     * Coefficient of the \f$1/\epsilon^2\f$ singularity 
     */
    int eps2;

    /**
     * Coefficient of the \f$1/\epsilon^2\f$ singularity 
     */
    int eps1;

    /**
     *  Coefficient of the finite term
     */
    double finite;
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLODrellYanBase();
  //@}

public:

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

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
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
   *  Virtual functions for the NLO calculation
   */
  //@{
  /**
   *  Calculate of the full next-to-leading order weight
   */
  virtual double NLOWeight() const;

  /**
   * Virtual matrix element, to be implemented in the
   * inheriting classes. The method should return the
   * loop matrix element including it'sa singular terms
   * Assumes the matrix element has the form
   * \f[ \frac{\alpha_S}{2\pi}\frac{C_F}{\Gamma(1-\epsilon)}
   *     \left(\frac{4\pi\mu^2}{\hat s}\right)^epsilon
   *    \left(\frac{A}{\epsilon^2}+\frac{B}{\epsilon}+C}
   * \f]
   * and A, B and C are returned.
   */
  virtual Singular virtualME() const = 0;

  /**
   * The leading-order matrix element, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  virtual double loME(const cPDVector & particles,
		      const vector<Lorentz5Momentum> & momenta,
		      bool first=false) const = 0;

  /**
   * The real matrix element divided by \f$2 g_S^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  virtual double realME(const cPDVector & particles,
			const vector<Lorentz5Momentum> & momenta) const = 0;
  //@}

  /**
   *  Subtracted virtual contribution
   */
  virtual double subtractedVirtual() const;

  /**
   *  Subtracted real contribution
   */
  virtual double subtractedReal(pair<double,double> x, double z,
				double zJac,
				double oldqPDF, double newqPDF, double newgPDF,
				bool order) const;

  /**
   *  Kinematic variables for the real radiation
   */
  //@{
  /**
   *  First  variable
   */
  double zTilde() const {return zTilde_;}

  /**
   *  Second variable
   */
  double vTilde() const {return vTilde_;}

  /**
   *  Azimuthal angle
   */
  double phi() const {return phi_;}
  //@}

  /**
   *  Calculate of the collinear counterterms
   */
  //@{
  /**
   *  Quark collinear counter term
   */
  double collinearQuark(double x, Energy2 mu2, double jac, double z,
			double oldPDF, double newPDF) const;

  /**
   *  Gluon collinear counter term
   */
  double collinearGluon(Energy2 mu2, double jac, double z,
			double oldPDF, double newPDF) const;
  //@}

  /**
   *   Dipole subtracted matrix elements
   */
  //@{
  /**
   *  \f$q\bar q\f$
   */
  double subtractedMEqqbar(const vector<Lorentz5Momentum> & pnew, bool order) const;

  /**
   *  \f$g\bar q\f$
   */
  double subtractedMEgqbar(const vector<Lorentz5Momentum> & pnew, bool order) const;

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
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<NLODrellYanBase> initNLODrellYanBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NLODrellYanBase & operator=(const NLODrellYanBase &);

private:

  /**
   *  Kinematic variables for the real radiation
   */
  //@{
  /**
   *  First  variable
   */
  mutable double zTilde_;

  /**
   *  Second variable
   */
  mutable double vTilde_;

  /**
   *  Azimuthal angle
   */
  mutable double phi_;
  //@}

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Pointer to the gluon ParticleData object
   */
  tcPDPtr gluon_;

  /**
   *  Factor for \f$C_F\f$ dependent pieces
   */
  mutable double CFfact_;

  /**
   *  Factor for \f$T_R\f$ dependent pieces
   */
  mutable double TRfact_;

  /**
   *  Strong coupling
   */
  mutable double alphaS_;

  /**
   *  Use a fixed value of \f$\alpha_S\f$
   */
  bool fixedAlphaS_;

  /**
   *  Leading-order matrix element
   */
  mutable double loME_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NLODrellYanBase. */
template <>
struct BaseClassTrait<Herwig::NLODrellYanBase,1> {
  /** Typedef of the first base class of NLODrellYanBase. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NLODrellYanBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NLODrellYanBase>
  : public ClassTraitsBase<Herwig::NLODrellYanBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NLODrellYanBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NLODrellYanBase is implemented. It may also include several, space-separated,
   * libraries if the class NLODrellYanBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwSusyNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NLODrellYanBase_H */
