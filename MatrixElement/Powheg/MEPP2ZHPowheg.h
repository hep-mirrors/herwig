// -*- C++ -*-
#ifndef HERWIG_MEPP2ZHPowheg_H
#define HERWIG_MEPP2ZHPowheg_H
//
// This is the declaration of the MEPP2ZHPowheg class.
//

#include "Herwig/MatrixElement/Hadron/MEPP2ZH.h"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2ZHPowheg class implements the matrix element
 * for \f$q\bar{q}\to Z^0h^0\f$.
 *
 * @see \ref MEPP2ZHPowhegInterfaces "The interfaces"
 * defined for MEPP2ZHPowheg.
 */
class MEPP2ZHPowheg: public MEPP2ZH {

public:

  /**
   * The default constructor.
   */
  MEPP2ZHPowheg();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
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
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   * Calculate the variable \f$x=M_{B}^2/s\f$ from the integration variables. 
   */
  double x(double xt, double v) const;

  /**
   * Calculate the momentum fraction of the first parton. 
   */
  double x_a(double x, double v) const;

  /**
   * Calculate the momentum fraction of second parton. 
   */
  double x_b(double x, double v) const;

  /**
   * Calculate the minimum of \f$x\f$. 
   */
  double xbar(double v) const;

  /**
   * Calculate the ratio of the radiative luminosity funcion to the
   * Born luminosity function for the \f$qg\f$ initiated channel. 
   */
  double Ltilde_qg(double x, double v) const;

  /**
   * Calculate the ratio of the radiative luminosity funcion to the
   * Born luminosity function for the \f$g\bar{q}\f$ initiated channel. 
   */
  double Ltilde_gq(double x, double v) const;

  /**
   * Calculate the ratio of the radiative luminosity funcion to the
   * Born luminosity function for the \f$q\bar{q}\f$ initiated channel. 
   */
  double Ltilde_qq(double x, double v) const;

  /**
   * Calculate the soft-virtual contribution to the NLO weight. 
   */
  double Vtilde_qq() const;

  /**
   * Function for calculation of the \f$g\bar{q}\f$ and \f$g\bar{q}\f$ 
   * initiated real contribution.
   */
  double Ccalbar_qg(double x) const;

  /**
   * Function for calculation of the \f$qg\f$ 
   * initiated real contribution.
   */
  double Fcal_qg(double x, double v) const;

  /**
   * Function for calculation of the \f$g\bar{q}\f$ initiated real
   * contribution.
   */
  double Fcal_gq(double x, double v) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Fcal_qq(double x, double v) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Ftilde_qg(double xt, double v) const;

  /**
   * Function for calculation of the \f$g\bar{q}\f$ initiated real
   * contribution.
   */
  double Ftilde_gq(double xt, double v) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Ftilde_qq(double xt, double v) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Ctilde_qg(double x, double v) const;

  /**
   * Function for calculation of the \f$g\bar{q}\f$ initiated real
   * contribution.
   */
  double Ctilde_gq(double x, double v) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Ctilde_qq(double x, double v) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
  static ClassDescription<MEPP2ZHPowheg> initMEPP2ZHPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2ZHPowheg & operator=(const MEPP2ZHPowheg &);

private:

  /**
   *  The momentum fraction of the first parton in the Born process
   */
  mutable double _xb_a;

  /**
   *  The momentum fraction of the second parton in the Born process
   */
  mutable double _xb_b;

  /**
   *  The ParticleData object for the first parton in the Born process
   */
  mutable tcPDPtr _parton_a;

  /**
   *  The ParticleData object for the second parton in the Born process
   */
  mutable tcPDPtr _parton_b;

  /**
   *  The BeamParticleData object for the first  hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer _hadron_A;

  /**
   *  The BeamParticleData object for the second hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer _hadron_B;

  /**
   *  the ParticleData object for the gluon
   */
  tcPDPtr _gluon;

  /**
   * The \f$T_R\f$ colour factor
   */
  const double TR_;

  /**
   *  The \f$C_F\f$ colour factor
   */
  const double CF_;

  /**
   *  The value of \f$\frac{\alpha_S}{2\pi}\f$ used for the calculation
   */
  mutable double _alphaS2Pi;

  /**
   *  The mass squared of the lepton pair
   */
  mutable Energy2 _mll2;

  /**  
   * The renormalization/factorization scale
   */
  mutable Energy2 _mu2;

  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int _contrib;

  /**
   *  Whether to use a fixed or a running QCD coupling for the NLO weight
   */
  unsigned int _nlo_alphaS_opt;

  /**
   *  The value of alphaS to use for the nlo weight if _nloalphaSopt=1
   */
  double _fixed_alphaS;

  /**
   *  The magnitude of the correction term to reduce the negative contribution
   */
  double _a;

  /**
   *  The power of the correction term to reduce the negative contribution
   */
  double _p;

  /**
   *  Cut-off for the correction function
   */
  double _eps;
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

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double _xt;

  /**
   *  The \f$v\f$ angular variable
   */
  double _v;
  //@}

  /**
   *  Values of the PDF's before radiation
   */
  //@{
  /**
   *  For the quark
   */
  mutable double _oldq;

  /**
   *  For the antiquark
   */
  mutable double _oldqbar;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2ZHPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2ZHPowheg,1> {
  /** Typedef of the first base class of MEPP2ZHPowheg. */
  typedef Herwig::MEPP2ZH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2ZHPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2ZHPowheg>
  : public ClassTraitsBase<Herwig::MEPP2ZHPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2ZHPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2ZHPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2ZHPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2ZHPowheg_H */
