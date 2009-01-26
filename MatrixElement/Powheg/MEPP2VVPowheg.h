// -*- C++ -*-
#ifndef HERWIG_MEPP2VVPowheg_H
#define HERWIG_MEPP2VVPowheg_H
//
// This is the declaration of the MEPP2VVPowheg class.
//

#include "Herwig++/MatrixElement/Hadron/MEPP2VV.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig++/MatrixElement/Powheg/NLO2to2Kinematics.h"
#include "Herwig++/Utilities/Maths.h"

namespace Herwig {
using namespace ThePEG;
using Math::ReLi2;
using Constants::pi;

/**
 * Here is the documentation of the MEPP2VVPowheg class.
 *
 * @see \ref MEPP2VVPowhegInterfaces "The interfaces"
 * defined for MEPP2VVPowheg.
 */
class MEPP2VVPowheg: public MEPP2VV {

public:

  /**
   * The default constructor.
   */
  MEPP2VVPowheg();

public:

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

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * This member check the collinear limits of the 
   * real emission matrix elements are equal to the 
   * appropriate combinations of Born ME's multiplied
   * by the splitting functions.
   */
  void sanityCheck() const;

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

  /**
   * Function to set the born variables. 
   */
  void getKinematics(double xt, double y);

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   * Calculate the ratio of the NLO luminosity to the LO
   * luminosity function for the \f$q\bar{q}\f$ initiated channel. 
   */
  double Lhat_ab(tcPDPtr a, tcPDPtr b, real2to3Kinematics Kinematics) const;

  /**
   * Calculate the universal soft-virtual contribution to the NLO weight. 
   */
  double Vtilde_universal(real2to3Kinematics S) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_qq_on_x(tcPDPtr a,tcPDPtr b,real2to3Kinematics C) const;

  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,real2to3Kinematics C) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qqb_on_x(tcPDPtr a,tcPDPtr b) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qg_on_x(tcPDPtr a,tcPDPtr b) const;

  /**
   * Function for calculation of the \f$gqb\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_gqb_on_x(tcPDPtr a,tcPDPtr b) const;

  /**
   * The regular part of the virtual correction matrix element(s).
   * For WZ production this is given by Equation B.2 in NPB 383 (1992) 
   * 3-44 *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_V_regular(real2to3Kinematics S) const;

  /**
   * The matrix element q + qb -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqb(real2to3Kinematics R) const;

  /**
   * The matrix element q + g  -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg(real2to3Kinematics R) const;

  /**
   * The matrix element g + qb -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gqb(real2to3Kinematics R) const;

  /**
   * The Born matrix element as given in Equation 3.1 - 3.3 in NPB 383 
   * (1992) *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_Born(born2to2Kinematics B) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2VVPowheg> initMEPP2VVPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VVPowheg & operator=(const MEPP2VVPowheg &);

private:

  /**
   *  Parameters for the NLO weight
   */
  //@{

  /**
   *  The BeamParticleData object for the plus and minus direction hadrons
   */
  Ptr<BeamParticleData>::transient_const_pointer hadron_A_;
  Ptr<BeamParticleData>::transient_const_pointer hadron_B_;

  /**
   * Born / virtual 2->2 kinematics.
   */
  born2to2Kinematics B_;

  /**
   * Soft limit of the 2->3 real emission kinematics.
   */
  real2to3Kinematics S_;

  /**
   * Soft-collinear limit of the 2->3 kinematics (emission in +z direction).
   */
  real2to3Kinematics SCp_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in -z direction).
   */
  real2to3Kinematics SCm_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in +z direction).
   */
  real2to3Kinematics Cp_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in -z direction).
   */
  real2to3Kinematics Cm_;

  /**
   * The resolved 2->3 real emission kinematics:
   */
  real2to3Kinematics H_;

  /**
   *  The ParticleData object for the plus and minus lo partons
   */
  tcPDPtr ab_, bb_;

  /**
   *  Flag indicating if the q & qbar are flipped or not i.e. this
   *  is true if q enters from the -z direction in the lab frame.
   */
  bool flipped_;

  /**
   *  Values of the PDF's before radiation
   */
  double lo_lumi_;

  /**
   *  The value of the leading order qqbar->VV matrix element
   */
  mutable double lo_me2_;

  /**
   *  The CF_ colour factor
   */
  double CF_;

  /**
   *  The TR_ colour factor
   */
  double TR_;

  /**
   *  The number of colours
   */
  double NC_;

  /**
   *  The up and down, left handed, quark-boson couplings
   */
  mutable double guL_, gdL_;

  /**
   *  The TGC coupling
   */
  mutable double eZ_;

  /**
   *  The CKM factor (Fij^2)
   */
  mutable double Fij2_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Whether to use a fixed or a running QCD coupling for the NLO weight
   */
  unsigned int nlo_alphaS_opt_;

  /**
   *  The value of alphaS to use for the nlo weight if nlo_alphaS_opt_=1
   */
  double fixed_alphaS_;

  /**
   *  Flag to remove or multiply in MCFM branching fractions for testing
   */
  unsigned int removebr_;
  //@}

  /**
   *  The value of \f$\alpha_S\f$ used for the calculation
   */
  mutable double alphaS_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VVPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2VVPowheg,1> {
  /** Typedef of the first base class of MEPP2VVPowheg. */
  typedef Herwig::MEPP2VV NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VVPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VVPowheg>
  : public ClassTraitsBase<Herwig::MEPP2VVPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VVPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VVPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VVPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegME.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VVPowheg_H */
