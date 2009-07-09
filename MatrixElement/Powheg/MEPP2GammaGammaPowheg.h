// -*- C++ -*-
#ifndef HERWIG_MEPP2GammaGammaPowheg_H
#define HERWIG_MEPP2GammaGammaPowheg_H
//
// This is the declaration of the MEPP2GammaGammaPowheg class.
//

#include "Herwig++/MatrixElement/Hadron/MEPP2GammaGamma.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2GammaGammaPowheg class.
 *
 * @see \ref MEPP2GammaGammaPowhegInterfaces "The interfaces"
 * defined for MEPP2GammaGammaPowheg.
 */
class MEPP2GammaGammaPowheg: public Herwig::MEPP2GammaGamma {

public:

  /**
   * The default constructor.
   */
  MEPP2GammaGammaPowheg();

public:

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
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
  static ClassDescription<MEPP2GammaGammaPowheg> initMEPP2GammaGammaPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2GammaGammaPowheg & operator=(const MEPP2GammaGammaPowheg &);

protected:

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   *  Dipole subtracted real matrix element for \f$q\bar q \to \gamma\gamma g\f$
   */
  double MEqqbarg(const vector<Lorentz5Momentum> &, bool) const;

  /**
   *  Dipole subtracted real matrix element for \f$q g \to \gamma\gamma q\f$
   */
  double MEqgq(const vector<Lorentz5Momentum> &) const;

  /**
   *  Dipole subtracted real matrix element for \f$g\bar q \to \gamma\gamma \bar q\f$
   */
  double MEqbargqbar(const vector<Lorentz5Momentum> &) const;

private:

  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;
  //@}

  /**
   *  Choice of the scale
   */
  //@{
  /**
   *  Type of scale
   */
  unsigned int scaleopt_;

  /**
   *  Fixed scale if used
   */
  Energy fixedScale_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;
  //@}

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

//   /**
//    *  Jacobian for \f$x_p\f$ integral
//    */
//   pair<double,double> jac_;

  /**
   *  Radiative variables
   */
  //@{
  /**
   *  \f$\tilde{z}\f$
   */
  double ztilde_;

  /**
   *   \f$\tilde{v}\f$
   */
  double vtilde_;

  /**
   *  \f$\phi\f$
   */
  double phi_;
  //@}

  /**
   *  Pointer to the gluon ParticleData object
   */
  tcPDPtr gluon_;

  /**
   *  Pointers to the vertices
   */
  //@{
  /**
   *  The photon vertex
   */
  AbstractFFVVertexPtr QEDVertex_;

  /**
   *  The gluon vertex
   */
  AbstractFFVVertexPtr QCDVertex_;
  //@}

  /**
   *  Strong coupling
   */
  mutable double alphaS_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2GammaGammaPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2GammaGammaPowheg,1> {
  /** Typedef of the first base class of MEPP2GammaGammaPowheg. */
  typedef Herwig::MEPP2GammaGamma NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2GammaGammaPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2GammaGammaPowheg>
  : public ClassTraitsBase<Herwig::MEPP2GammaGammaPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2GammaGammaPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2GammaGammaPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2GammaGammaPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2GammaGammaPowheg_H */
