// -*- C++ -*-
#ifndef HERWIG_MEPP2VVPowheg_H
#define HERWIG_MEPP2VVPowheg_H
//
// This is the declaration of the MEPP2VVPowheg class.
//

#include "Herwig/MatrixElement/Hadron/MEPP2VV.h"
#include "Herwig/MatrixElement/Powheg/VVKinematics.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

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
 
  /** @name Member functions for the generation of hard QCD radiation */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return ISR;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr,
						 ShowerInteraction inter);
  //@}

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
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * This member check the collinear limits of the 
   * real emission matrix elements are equal to the 
   * appropriate combinations of Born ME's multiplied
   * by the splitting functions.
   */
  bool sanityCheck() const;

  /**
   * Return the CKM matrix elements.
   */
  Complex CKM(int ix,int iy) const { return ckm_[ix][iy]; }

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

public:

  /**
   * Function to set the born variables. 
   */
  void getKinematics(double xt, double y, double theta2);

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   * Calculate the ratio of the NLO luminosity to the LO
   * luminosity function for the \f$q\bar{q}\f$ initiated channel. 
   */
  double Lhat_ab(tcPDPtr a, tcPDPtr b, realVVKinematics Kinematics) const;

  /**
   * Calculate the universal soft-virtual contribution to the NLO weight. 
   */
  double Vtilde_universal(realVVKinematics S) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_qq_on_x(tcPDPtr a,tcPDPtr b,realVVKinematics C) const;

  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,realVVKinematics C) const;

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
  double M_V_regular(realVVKinematics S) const;

  /**
   *  Member variable to store the
   * regular part of the virtual correction matrix element(s).
   * For WZ production this is given by Equation B.2 in NPB 383 (1992) 
   * 3-44 *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element.
   */
  mutable double M_V_regular_;

  /**
   * The matrix element q + qb -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqb(realVVKinematics R) const;

  /**
   *  Member variable to store the matrix element q + qb -> n + g times tk*uk 
   */
  mutable Energy2 t_u_M_R_qqb_;

  /**
   * The matrix element q + g  -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg(realVVKinematics R) const;

  /**
   *  Member variable to store the matrix element q + g  -> n + q times tk*uk 
   */
  mutable Energy2 t_u_M_R_qg_;

  /**
   * The matrix element g + qb -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gqb(realVVKinematics R) const;

  /**
   *  Member variable to store the matrix element g + qb -> n + q times tk*uk 
   */
  mutable Energy2 t_u_M_R_gqb_;

  /**
   * The matrix element q + qb -> n + g times (tk*uk)^2 - using helicity amplitudes
   */
  Energy2 t_u_M_R_qqb_hel_amp(realVVKinematics R) const;

  /**
   *  Member variable to store the
   * matrix element q + qb -> n + g times (tk*uk)^2 - using helicity amplitudes
   */
  mutable Energy2 t_u_M_R_qqb_hel_amp_;

  /**
   * The matrix element q + g -> n + q times (tk*uk)^2 - using helicity amplitudes
   */
  Energy2 t_u_M_R_qg_hel_amp(realVVKinematics R) const;

  /**
   *  Member variable to store the
   * matrix element q + g -> n + q times (tk*uk)^2 - using helicity amplitudes
   */
  mutable Energy2 t_u_M_R_qg_hel_amp_;

  /**
   * The matrix element g + qb -> n + qb times (tk*uk)^2 - using helicity amplitudes
   */
  Energy2 t_u_M_R_gqb_hel_amp(realVVKinematics R) const;

  /**
   *  Member variable to store the
   * matrix element g + qb -> n + qb times (tk*uk)^2 - using helicity amplitudes
   */
  mutable Energy2 t_u_M_R_gqb_hel_amp_;

  /**
   * The leading order matrix element - using helicity amplitudes
   */
  double lo_me() const;

  /**
   * The Born matrix element as given in Equation 3.1 - 3.3 in NPB 383 
   * (1992) *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_Born_WZ(bornVVKinematics B) const;

  /**
   *  Member variable to store the 
   * Born matrix element as given in Equation 3.1 - 3.3 in NPB 383 
   * (1992) *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  mutable double M_Born_;

  /**
   * The Born matrix element as given in Equation 2.18 - 2.19 in NPB 357 
   * (1991) *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_Born_ZZ(bornVVKinematics B) const;

  /**
   * M_V_regular_ZZ is the regular part of the one-loop ZZ matrix element 
   * exactly as defined in Eqs. B.1 & B.2 of  NPB 357(1991)409-438 ***
   * modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_V_regular_ZZ(realVVKinematics S) const;

  /**
   * t_u_M_R_qqb_ZZ is the q + qb -> n + g times tk*uk real emission 
   * matrix element as defined in Eq. C.1 of  NPB 357(1991)409-438 ***
   * modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  Energy2 t_u_M_R_qqb_ZZ(realVVKinematics R) const;

  /**
   * The Born matrix element as given in Equation 3.2 - 3.8 in NPB 410 
   * (1993) *** modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_Born_WW(bornVVKinematics B) const;

  /**
   * M_V_regular_WW is the regular part of the one-loop WW matrix element 
   * exactly as defined in Eqs. C.1 - C.7 of of NPB 410(1993)280-324 ***
   * modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  double M_V_regular_WW(realVVKinematics S) const;

  /**
   * t_u_M_R_qqb_WW is the q + qb -> n + g times tk*uk real emission 
   * matrix element as defined in Eq. D.1-D.5 of  NPB 410(1993)280-324 ***
   * modulo a factor 1/(2s) ***, which is a flux factor that 
   * those authors absorb in the matrix element. 
   */
  Energy2 t_u_M_R_qqb_WW(realVVKinematics R) const;

  /**
   * Return the factorion scale squared.
   */
  Energy2 mu_F2() const;

  /**
   * Return the renormalisation scale squared.
   */
  Energy2 mu_UV2() const;

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
  MEPP2VVPowheg & operator=(const MEPP2VVPowheg &) = delete;

private:

  /**
   *  Parameters for the NLO weight
   */
  //@{

  /**
   * Parameter to determine when to use limiting value of real emission
   * matrix elements, to avoid rounding error issues.
   */
  double tiny;

  /**
   *  The BeamParticleData object for the plus and minus direction hadrons
   */
  tcBeamPtr hadron_A_;
  tcBeamPtr hadron_B_;

  /**
   * Born / virtual 2->2 kinematics.
   */
  bornVVKinematics B_;

  /**
   * Soft limit of the 2->3 real emission kinematics.
   */
  realVVKinematics S_;

  /**
   * Soft-collinear limit of the 2->3 kinematics (emission in +z direction).
   */
  realVVKinematics SCp_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in -z direction).
   */
  realVVKinematics SCm_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in +z direction).
   */
  realVVKinematics Cp_;

  /**
   * The collinear limit of the 2->3 kinematics (emission in -z direction).
   */
  realVVKinematics Cm_;

  /**
   * The resolved 2->3 real emission kinematics:
   */
  realVVKinematics H_;

  /**
   *  The ParticleData object for the plus and minus lo partons
   */
  tcPDPtr ab_, bb_;

  /**
   *  The ParticleData object for the quark and antiquark 
   *  (which can be in a different order to ab_ and bb_).
   */
  tcPDPtr quark_, antiquark_;

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
   *  The weak coupling and the sin (squared) of the Weinberg angle
   */
  mutable double gW_, sin2ThetaW_;

  /**
   *  The up and down, left handed, quark-boson couplings
   */
  mutable double guL_, gdL_;

  /**
   *  The up and down, right handed, quark-boson couplings (for WW & ZZ)
   */
  mutable double guR_, gdR_;

  /**
   *  The TGC coupling
   */
  mutable double eZ_;

  /**
   *  The TGC coupling squared. This is useful for debugging purposes
   *  when one wants to turn of t-channel * TGC interference contributions
   *  but leave the pure TGC contributions intact. It is also needed in 
   *  order to transform WZ matrix elements into WW ones.
   */
  mutable double eZ2_;

  /**
   *  The CKM factor (Fij^2)
   */
  mutable double Fij2_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Whether to generate all channels contributions or just qqb or just 
   *  qg+gqb contributions
   */
  unsigned int channels_;

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

  /**
   * Selects a dynamic (sHat) or fixed factorization scale
   */
  unsigned int scaleopt_;

  /**
   * The factorization scale
   */
  Energy mu_F_;

  /**
   * The renormalization scale
   */
  Energy mu_UV_;

  /**
   * The pT of V1 in a radiative event in the lab frame (for scale setting only)
   */
  Energy2 k1r_perp2_lab_;

  /**
   * The pT of V2 in a radiative event in the lab frame (for scale setting only)
   */
  Energy2 k2r_perp2_lab_;

  /**
   * The ckm matrix elements (unsquared, to allow interference)
   */
  vector< vector<Complex> > ckm_;

  /**
   * Option to impose helicity conservation on the real NLO ME's (greatly improves evaluation time).
   */
  bool helicityConservation_;

  /**
   *  The q + qb -> v1 + v2 + g  helicity amplitudes  
   */
  mutable ProductionMatrixElement qqb_hel_amps_;

  /**
   *  The q + g  -> v1 + v2 + q  helicity amplitudes  
   */
  mutable ProductionMatrixElement qg_hel_amps_;

  /**
   *  The g + qb -> v1 + v2 + qb helicity amplitudes  
   */
  mutable ProductionMatrixElement gqb_hel_amps_;
  //@}

  /**
   *  The vertices
   */
  //@{
  /**
   *  The photon fermion-antifermion vertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   *  The W fermion-antifermion vertex
   */
  AbstractFFVVertexPtr FFWvertex_;

  /**
   *  The Z fermion-antifermionvertex
   */
  AbstractFFVVertexPtr FFZvertex_;

  /**
   *  The triple electroweak gauge boson vertex
   */
  AbstractVVVVertexPtr WWWvertex_;

  /**
   *  The quark-antiquark gluon vertex
   */
  AbstractFFVVertexPtr FFGvertex_;
  //@}

  /**
   *  The value of \f$\alpha_S\f$ used for the calculation
   */
  mutable double alphaS_;

protected:

  /**
   * Returns the matrix element for a given type of process,
   * rapidity of the jet \f$y_j\f$ and transverse momentum \f$p_T\f$
   * @param emis_type the type of emission,
   * (0 is \f$q\bar{q}\to Vg\f$, 1 is \f$qg\to Vq\f$ and 2 is \f$g\bar{q}\to V\bar{q}\f$)
   * @param pT The transverse momentum of the jet
   * @param R The object containing the kinematics
   */
  double getResult(int emis_type, realVVKinematics R, Energy pT);
 
  /**
   *  generates the hardest emission (yj,p)
   * @param pnew The momenta of the new particles
   * @param emissiontype The type of emission, as for getResult
   * @return Whether not an emission was generated
   */
  bool getEvent(vector<Lorentz5Momentum> & pnew,unsigned int & emissiontype);
  
  /**
   *  sets the QCD, EW and PDF scales
   * @param pT The pT of the current step in the veto algorithm
   */
  void setTheScales(Energy pT);

  /**
   * The matrix element q + qb -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqb_hel_amp(realVVKinematics R, bool getMatrix) const;


  /**
   * The matrix element q + g  -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg_hel_amp(realVVKinematics R, bool getMatrix) const;

  /**
   * The matrix element g + qb -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gqb_hel_amp(realVVKinematics R, bool getMatrix) const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  double lo_me(bool getMatrix) const;

  /**
   * Recalculate hard vertex to include spin correlations for radiative events.
   */
  void recalculateVertex();

  /**
   * Member which selects a two body decay mode for each vector
   * boson and distributes decay products isotropically
   */
  bool isotropicDecayer();

  /**
   * The triangle function lambda(x,y,z)=sqrt(x^2+y^2+z^2-2*x*y-2*y*z-2*x*z)
   */
  Energy2 triangleFn(Energy2,Energy2,Energy2);

private:

  /**
   * If this boolean is true the n+1 body helicity amplitudes will be
   * used to calculate a hard vertex based on those kinematics for spin
   * correlations in the decays.
   */
  bool realMESpinCorrelations_;

  /**
   * The colour & spin averaged n-body (leading order) matrix element squared.
   */
  double lo_me_;

  /**
   * The resolved 2->3 real emission kinematics.
   */
  realVVKinematics R_;

  /**
   * This specifies the emitting configuration: 
   * 1: q + qbar -> V1 + V2 + g
   * 2: q + g    -> V1 + V2 + q
   * 3: g + qbar -> V1 + V2 + qbar.
   */
  unsigned int channel_;

  /**
   * Identifies the space-like mother of the branching
   * as quark (+1) or antiquark (-1):
   */
  int fermionNumberOfMother_;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr showerAlphaS_;

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double power_;

  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preqqbar_;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double preqg_;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double pregqbar_;

  /**
   * The QCD beta function divided by 4pi, (11-2/3*nf)/4/pi, with nf = 5.
   */
  double b0_;

  /**
   * The fundamental QCD scale in the one-loop alpha_{S} used for the crude
   * (not the very crude) overestimate of the Sudakov exponent. The default
   * value is set so such that alphaS(MZ), neglecting all flavour threshold
   * effects i.e. MZ*exp(-1/2/b0_/alphaS(MZ)).
   */
  Energy LambdaQCD_;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> prefactor_;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the ShowerProgenitor objects for the partons
   */
  PPtr qProgenitor_;
  PPtr qbProgenitor_;

  /**
   *  Pointers to the Shower particle objects for the partons
   */
  PPtr showerQuark_;
  PPtr showerAntiquark_;

  /**
   *  Pointers to the BeamParticleData objects
   */
  tcBeamPtr qHadron_;
  tcBeamPtr qbHadron_;
  //@}

  /**
   *  Properties of the boson and jets
   */
  //@{
  /**
   *  Pointers to the Shower particle objects for the partons
   */
  PPtr gluon_;
  PPtr V1_;
  PPtr V2_;
  vector<PPtr> children_;

  /**
   *  Flag indicating if the q & qbar are flipped or not i.e. this
   *  is true if q enters from the -z direction in the lab frame.
   */
  bool flipped_;

  /**
   *  the rapidity of the jet
   */
  double Yk_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pT_;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy min_pT_;

  // Work out the scales we want to use in the matrix elements and the pdfs:
  /**
   * Scale for alpha_S: pT^2 of the diboson system.
   */
  Energy2 QCDScale_;

  /**
   * Scale for real emission PDF: 
   */
  Energy2 PDFScale_;

  /**
   * Scale of electroweak vertices: mVV^2 the invariant mass of the diboson system.
   */
  Energy2 EWScale_;

  /**
   * A matrix to hold the home-grown production matrix element
   */
  mutable Complex productionMatrix_[3][3][3][3];

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
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VVPowheg_H */
