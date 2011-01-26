// -*- C++ -*-
#ifndef HERWIG_VVKinematics_H
#define HERWIG_VVKinematics_H
//
// This is the declaration of the VVKinematics class.
//

#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup VVKinematics
 *  The bornVVKinematics class is used to store information on the 
 *  the kinematics of the real emission processes needed for the 
 *  evaluation of matrix elements in the real part of the NLO process.
 */
class bornVVKinematics {

public:

  /**
   * Default constructor 
   */
  bornVVKinematics();

  /**
   * Meaningful constructor: takes the momenta 
   * and Bjorken x values of the partons involved
   * in the 2->2 scattering of the leading order
   * / virtual process and calculates all interesting
   * Mandelstam and Born variables.
   */
  bornVVKinematics(vector<Lorentz5Momentum> Momenta, double x1, double x2);

public:

  /**
   * Read-only access to all of the above member variables.
   */

  /**
   * @name Leading order momentum fractions and associated etabar's:
   */
  //@{
  /**
   * Leading order momentum fraction for first particle
   */
  double x1b()   const { return x1b_; }

  /**
   * Leading order etabar for first particle
   */
  double eta1b() const { return eta1b_; } 

  /**
   * Leading order momentum fraction for second particle
   */
  double x2b()   const { return x2b_; }

  /**
   * Leading order etabar for second particle
   */
  double eta2b() const { return eta2b_; }
  //@}

  /**
   * @name The Born momenta according to the notation of the FMNR papers,
   * in the diboson centre of mass frame:
   */
  //@{
  /**
   * Momentum \f$p_1\f$
   */
  Lorentz5Momentum p1b() const { return p1b_; }
  
  /**
   * Momentum \f$p_2\f$
   */
  Lorentz5Momentum p2b() const { return p2b_; }
  
  /**
   * Momentum \f$k_1\f$
   */
  Lorentz5Momentum k1b() const { return k1b_; }
  
  /**
   * Momentum \f$k_2\f$
   */
  Lorentz5Momentum k2b() const { return k2b_; }
  //@}

  /**
   * @name The diboson invariant mass / shat, that and uhat:
   */
  //@{
  /**
   * \f$\hat s\f$
   */
  Energy2 sb() const { return sb_; }

  /**
   * \f$\hat t\f$
   */
  Energy2 tb() const { return tb_; }

  /**
   * \f$\hat u\f$
   */
  Energy2 ub() const { return ub_; }
  //@}

  /**
   * The diboson rapidity:
   * Note Yb_ = + lastY() if flipped = false 
   * but  Yb_ = - lastY() if flipped = true.
   * Yb_ is always defined with the quark travelling in the +z direction!
   */
  double Yb() const { return Yb_; }

  /**
   * @name Masses of the final state bosons:
   */
  //@{
  /**
   *  Mass squared og the first boson
   */
  Energy2 k12b() const { return k12b_; }

  /**
   *  Mass squared og the second boson
   */
  Energy2 k22b() const { return k22b_; }
  //@}

  /**
   * Polar and azimuthal angles of the dibosons in their rest frame:
   */
  double theta1b() const { return theta1b_; }

  /**
   * A check to make sure that the momenta calculated from 
   * the energies and angles are equal to those of meMomenta().
   */
  void sanityCheck() const;

 private:

  /**
   * Invariants required for the evaluation of 2-> 2 next-to-leading 
   * order quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */

  /**
   * @name Leading order momentum fractions and associated etabar's:
   */
  //@{
  /**
   * Leading order momentum fraction for first particle
   */
  double x1b_;

  /**
   * Leading order etabar for first particle
   */
  double eta1b_;

  /**
   * Leading order momentum fraction for second particle
   */
  double x2b_;

  /**
   * Leading order etabar for second particle
   */
  double eta2b_;
  //@}

  /**
   * @name The Born momenta according to the notation of the FMNR papers,
   * in the diboson centre of mass frame:
   */
  //@{
  /**
   * Momentum \f$p_1\f$
   */
  Lorentz5Momentum p1b_;
  
  /**
   * Momentum \f$p_2\f$
   */
  Lorentz5Momentum p2b_;
  
  /**
   * Momentum \f$k_1\f$
   */
  Lorentz5Momentum k1b_;
  
  /**
   * Momentum \f$k_2\f$
   */
  Lorentz5Momentum k2b_;
  //@}

  /**
   * @name The diboson invariant mass / shat, that and uhat:
   */
  //@{
  /**
   * \f$\hat s\f$
   */
  Energy2 sb_;

  /**
   * \f$\hat t\f$
   */
  Energy2 tb_;

  /**
   * \f$\hat u\f$
   */
  Energy2 ub_;
  //@}

  /**
   * The diboson rapidity:
   * Note Yb_ = + lastY() if flipped = false 
   * but  Yb_ = - lastY() if flipped = true.
   * Yb_ is always defined with the quark travelling in the +z direction!
   */
  double Yb_;

  /**
   * @name Masses of the final state bosons:
   */
  //@{
  /**
   *  Mass squared og the first boson
   */
  Energy2 k12b_;

  /**
   *  Mass squared og the second boson
   */
  Energy2 k22b_;
  //@}

  /**
   * Polar angle of the dibosons in their rest frame:
   */
  double theta1b_;     
};



/** \ingroup VVKinematics
 *  The realVVKinematics class is used to store information on the 
 *  the kinematics of the real emission processes needed for the 
 *  evaluation of matrix elements in the real part of the NLO process.
 */
class realVVKinematics {

public:

  /**
   * Default constructor 
   */
  realVVKinematics();

  /**
   * Meaningful constructor: takes the Born variables
   * from the leading order /virtual 2->2 process and 
   * the raw \f$\tilde{x}, y\f$ radiative variables and turns
   * these into a set of 2->3 momenta with associated 
   * Mandelstam variables, Bjorken x values etc. 
   * @param bornVariables The object for the Born kinematics
   * @param xt The \f$\tilde{x}\f$ radiative variable.
   * @param y  The angular radiative variable (the cosine 
   * of the  polar angle of the emitted gluon in the partonic  
   * @param theta2 The angle \f$\theta_2\f$
   * CMS frame). 
   */
  realVVKinematics(bornVVKinematics bornVariables,double xt, double y, double theta2);

  /**
   * A check to make sure that the momenta calculated from 
   * the energies and angles are equal to those of meMomenta().
   */
  void sanityCheck() const;
     
public:

  /**
   * Read-only access to all of the above member variables.
   */

  /**
   * The bornVVKinematics underlying the 2->3 kinematics
   */
  bornVVKinematics bornVariables() const { return bornVariables_; }

  /**
   * The lower bound on the x integration:
   */
  double xbar() const { return xbar_; }

  /**
   * @name The `raw' radiative variables.
   */
  //@{
  /**
   *  The \f$\tilde{x}\f$ radiative variable
   */
  double xt() const { return xt_; }

  /**
   *  The \f$y\f$ radiative variable
   */
  double y() const { return y_; }

  /**
   *  The \f$x_r\f$ radiative variable
   */
  double xr() const { return xr_; }
  //@}

  /**
   * The momentum fraction of the parton incident from the +z direction.
   */
  double x1r() const { return x1r_; }

  /**
   * The momentum fraction of the parton incident from the -z direction.
   */
  double x2r() const { return x2r_; }

  /**
   * Invariants required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */

  /**
   * @name First the Born variables:
   */
  //@{
  /**
   *  \f$s_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 s2r() const { return s2r_; }

  /**
   *  \f$k_1^2\f$ mass of first vector boson
   */
  Energy2 k12r() const { return k12r_; }

  /**
   *  \f$k_2^2\f$ mass of second vector boson
   */
  Energy2 k22r() const { return k22r_; }

  /**
   * \f$\theta_1\f$ angle from Frixione et al. NPB.383,3
   */
  double  theta1r() const { return theta1r_; }

  /**
   * \f$\theta_2\f$ angle from Frixione et al. NPB.383,3
   */
  double  theta2r() const { return theta2r_; }
  //@}

  /**
   * @name Then the rest:
   */
  //@{
  /**
   *  \f$p_T^2\f$ in the lab frame
   */
  Energy2 pT2_in_lab() const { return tkr_*ukr_/sr_; }

  /**
   * \f$s\f$ from Frixione et al. NPB.383,3
   */
  Energy2 sr() const { return sr_; }

  /**
   * \f$t_k\f$ from Frixione et al. NPB.383,3
   */
  Energy2 tkr() const { return tkr_; }   

  /**
   * \f$u_k\f$ from Frixione et al. NPB.383,3
   */
  Energy2 ukr() const { return ukr_; }

  /**
   * \f$\cos\psi\f$ from Frixione et al. NPB.383,3
   */
  double  cpsir() const { return cpsir_; }

  /**
   * \f$\cos\psi'\f$ from Frixione et al. NPB.383,3
   */
  double  cpsipr() const { return cpsiprr_; }

  /**
   * \f$\beta_x\f$ from Frixione et al. NPB.383,3
   */
  double  betaxr() const { return betaxr_; }

  /**
   * \f$v_1\f$ variable from Frixione et al. NPB.383,3
   */
  double  v1r() const { return v1r_; }

  /**
   * \f$v_2\f$ variable from Frixione et al. NPB.383,3
   */
  double  v2r() const { return v2r_; }

  /**
   * \f$q_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q1r() const { return q1r_; }

  /**
   * \f$q_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q2r() const { return q2r_; } 

  /**
   * \f$\hat q_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q1hatr() const { return q1hatr_; }

  /**
   * \f$\hat q_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q2hatr() const { return q2hatr_; }

  /**
   * \f$\hat w_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 w1r() const { return w1r_; }

  /**
   * \f$\hat w_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 w2r() const { return w2r_; }

  /**
   *  4-momentum \f$p_1\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum p1r() const { return p1r_; }

  /**
   *  4-momentum \f$p_2\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum p2r() const { return p2r_; }

  /**
   *  4-momentum \f$k\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum kr()  const { return kr_ ; }

  /**
   *  4-momentum \f$k_1\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum k1r() const { return k1r_; }

  /**
   *  4-momentum \f$k_2\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum k2r() const { return k2r_; }

  /**
   *  Set 4-momentum \f$p_1\f$ from Frixione et al. NPB.383,3
   */
  void p1r(Lorentz5Momentum p1r) { p1r_ = p1r; }

  /**
   *  Set 4-momentum \f$p_2\f$ from Frixione et al. NPB.383,3
   */
  void p2r(Lorentz5Momentum p2r) { p2r_ = p2r; }

  /**
   *  Set 4-momentum \f$k\f$ from Frixione et al. NPB.383,3
   */
  void kr (Lorentz5Momentum kr ) { kr_  = kr ; }

  /**
   *  Set 4-momentum \f$k_1\f$ from Frixione et al. NPB.383,3
   */
  void k1r(Lorentz5Momentum k1r) { k1r_ = k1r; }

  /**
   *  Set 4-momentum \f$k_2\f$ from Frixione et al. NPB.383,3
   */
  void k2r(Lorentz5Momentum k2r) { k2r_ = k2r; }
  //@}
 
private:

  /**
   * The bornVVKinematics object underlying the 2->3 kinematics.
   */
  bornVVKinematics bornVariables_;

  /**
   * The lower bound on the x integration.
   */
  double xbar_;

  // The `raw' radiative variables.
  /**
   * @name The `raw' radiative variables.
   */
  //@{
  /**
   *  The \f$\tilde{x}\f$ radiative variable
   */
  double xt_;

  /**
   *  The \f$y\f$ radiative variable
   */
  double y_;

  /**
   *  The \f$x_r\f$ radiative variable
   */
  double xr_;
  //@}

  /**
   * The momentum fraction of the parton incident from the +z direction.
   */
  double x1r_;

  /**
   * The momentum fraction of the parton incident from the -z direction.
   */
  double x2r_;

  /**
   * Invariants required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */

  /**
   * @name First the Born variables:
   */
  //@{
  /**
   *  \f$s_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 s2r_;

  /**
   *  \f$k_1^2\f$ mass of first vector boson
   */
  Energy2 k12r_;

  /**
   *  \f$k_2^2\f$ mass of second vector boson
   */
  Energy2 k22r_;

  /**
   * \f$\theta_1\f$ angle from Frixione et al. NPB.383,3
   */
  double  theta1r_;

  /**
   * \f$\theta_2\f$ angle from Frixione et al. NPB.383,3
   */
  double  theta2r_;
  //@}

  /**
   * @name Then the rest:
   */
  //@{
  /**
   * \f$s\f$ from Frixione et al. NPB.383,3
   */
  Energy2 sr_;

  /**
   * \f$t_k\f$ from Frixione et al. NPB.383,3
   */
  Energy2 tkr_;

  /**
   * \f$u_k\f$ from Frixione et al. NPB.383,3
   */
  Energy2 ukr_;

  /**
   * \f$\cos\psi\f$ from Frixione et al. NPB.383,3
   */
  double  cpsir_;

  /**
   * \f$\cos\psi'\f$ from Frixione et al. NPB.383,3
   */
  double  cpsiprr_;

  /**
   * \f$\beta_x\f$ from Frixione et al. NPB.383,3
   */
  double  betaxr_;

  /**
   * \f$v_1\f$ variable from Frixione et al. NPB.383,3
   */
  double  v1r_;

  /**
   * \f$v_2\f$ variable from Frixione et al. NPB.383,3
   */
  double  v2r_;

  /**
   * \f$q_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q1r_;

  /**
   * \f$q_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q2r_;

  /**
   * \f$\hat q_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q1hatr_;

  /**
   * \f$\hat q_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 q2hatr_;

  /**
   * \f$\hat w_1\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 w1r_;

  /**
   * \f$\hat w_2\f$ variable from Frixione et al. NPB.383,3
   */
  Energy2 w2r_;

  /**
   *  4-momentum \f$p_1\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum p1r_;

  /**
   *  4-momentum \f$p_2\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum p2r_;

  /**
   *  4-momentum \f$k\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum kr_;

  /**
   *  4-momentum \f$k_1\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum k1r_;

  /**
   *  4-momentum \f$k_2\f$ from Frixione et al. NPB.383,3
   */
  Lorentz5Momentum k2r_;
  //@}

};

}

#endif /* HERWIG_VVKinematics_H */
