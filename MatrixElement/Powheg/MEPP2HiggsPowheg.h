// -*- C++ -*-
//
// MEPP2HiggsPowheg.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2HiggsPowheg_H
#define HERWIG_MEPP2HiggsPowheg_H
//
// This is the declaration of the MEPP2HiggsPowheg class.
//

#include "Herwig/MatrixElement/Hadron/MEPP2Higgs.h"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
 
/**
 * The MEPP2HiggsPowheg class implements the matrix element for the process
 * pp->Higgs with different Higgs shape prescriptions (see details in hep-ph/9505211)
 * and the NLL corrected Higgs width (see details in the FORTRAN HERWIG manual).
 *
 * @see \ref MEPP2HiggsPowhegInterfaces "The interfaces"
 * defined for MEPP2HiggsPowheg.
 */
class MEPP2HiggsPowheg: public MEPP2Higgs {

public:

  /**
   * The default constructor.
   */
  MEPP2HiggsPowheg();

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
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

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

  /**
   * Function to set the born variables. 
   */
  void get_born_variables() const;

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   * Invariant required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */
  Energy2 s(double xt, double y) const {
    return  p2_/x(xt,y);
  }

  /**
   * Invariant required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */
  Energy2 tk(double xt, double y) const {
    double  x_xt_y(x(xt,y));
    return -0.5*p2_/x_xt_y*(1.- x_xt_y)*(1.-y);
  }

  /**
   * Invariant required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */
  Energy2 uk(double xt, double y) const {
    double  x_xt_y(x(xt,y));
    return -0.5*p2_/x_xt_y*(1.- x_xt_y)*(1.+y);
  }

  /**
   * Calculate the minimum of \f$x\f$. 
   */
  double xbar(double y) const;

  /**
   * Calculate auxiliary function of \f$\bar{x}(y)\f$, \f$\bar{\eta}(y)\f$. 
   */
  double etabar(double y) const;

  /**
   * Calculate the variable \f$x=p^{2}/s\f$ from the integration variables. 
   */
  double x(double xt, double y) const {
    double x0(xbar(y));
    return x0+(1.-x0)*xt;
  }

  /**
   * Calculate the momentum fraction of the plus parton. 
   */
  double xp(double x, double y) const;

  /**
   * Calculate the momentum fraction of the minus parton. 
   */
  double xm(double x, double y) const;

  /**
   * Calculate the ratio of the NLO luminosity to the LO
   * luminosity function for the \f$q\bar{q}\f$ initiated channel. 
   */
  double Lhat_ab(tcPDPtr a, tcPDPtr b, double x, double y) const;

  /**
   * Calculate the universal soft-virtual contribution to the NLO weight. 
   */
  double Vtilde_universal() const;

  /**
   * Function for calculation of the \f$gg\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_qg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * The regular part of the virtual correction matrix element(s) 
   */
  double M_V_regular() const;

  /**
   * The matrix element q + qbar -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqbar(double xt, double y) const;

  /**
   * The matrix element qbar + q -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qbarq(double xt, double y) const;

  /**
   * The matrix element g + g    -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_gg(double xt, double y) const;

  /**
   * The matrix element q + g    -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg(double xt, double y) const;

  /**
   * The matrix element g + q    -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gq(double xt, double y) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qqbar_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$\bar{q}q\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qbarq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$qq\f$ 
   * initiated real contribution.
   */
  double Rtilde_Ltilde_gg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsPowheg & operator=(const MEPP2HiggsPowheg &) = delete;

private:

  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  The colour factors
   */
  const double CF_ , CA_ , TR_;

  /**
   * Number of light flavours (in the beta function beta0_)
   */
  const int nlf_;

  /**
   * (Proportional to) The beta function
   */
  const double beta0_;

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
  //@}

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double xt_;

  /**
   *  The \f$y\f$ angular variable
   */
  double y_;
  //@}

  /**
   *  Values of the PDF's before radiation
   */
  mutable double lo_lumi_;

  /**
   *  The value of the leading order gg->H matrix element
   */
  mutable double lo_ggME_;

  /**
   * The invariant mass of the lo final state. 
   */
  mutable Energy2 p2_      ;

  /**
   * The invariant mass of the lo final state. 
   */
  mutable Energy2 s2_      ;

  /**
   *  The momentum fraction of the plus parton in the Born process
   */
  mutable double xbp_;

  /**
   *  The momentum fraction of the minus parton in the Born process
   */
  mutable double  xbm_;

  /**
   *  The sqrt(1-xbp_) 
   */
  mutable double etabarp_;

  /**
   *  The sqrt(1-xbm_) 
   */
  mutable double etabarm_;

  /**
   *  The ParticleData object for the plus lo parton
   */
  mutable tcPDPtr a_lo_;

  /**
   *  The ParticleData object for the minus lo parton
   */
  mutable tcPDPtr b_lo_;

  /**
   *  The BeamParticleData object for the plus direction hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer hadron_A_;

  /**
   *  The BeamParticleData object for the  minus direction hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer hadron_B_;

  /**
   *  The value of \f$\alpha_S\f$ used for the calculation
   */
  mutable double alphaS_;

  /**
   * Selects a dynamic (sHat) or fixed factorization scale
   */
  unsigned int scaleopt_;

  /**
   * The factorization  scale 
   */
  Energy mu_F_;

  /**
   * The renormalization scale
   */
  Energy mu_UV_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;
};

}

#endif /* HERWIG_MEPP2HiggsPowheg_H */
