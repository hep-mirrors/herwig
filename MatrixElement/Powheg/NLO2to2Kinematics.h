// -*- C++ -*-
#ifndef HERWIG_NLO2to2Kinematics_H
#define HERWIG_NLO2to2Kinematics_H
//
// This is the declaration of the NLO2to2Kinematics class.
//

#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup NLO2to2Kinematics
 *  The born2to2Kinematics class is used to store information on the 
 *  the kinematics of the real emission processes needed for the 
 *  evaluation of matrix elements in the real part of the NLO process.
 */
class born2to2Kinematics {
 
private:
  /**
   * Invariants required for the evaluation of 2-> 2 next-to-leading 
   * order quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */

  // Leading order momentum fractions and associated etabar's:
  double xpb_;
  double etapb_;
  double xmb_;
  double etamb_;

  // The Born momenta according to the notation of the FMNR papers,
  // in the diboson centre of mass frame:
  Lorentz5Momentum p1b_;
  Lorentz5Momentum p2b_;
  Lorentz5Momentum k1b_;
  Lorentz5Momentum k2b_;

  // The diboson invariant mass / shat, that and uhat:
  Energy2 sb_;
  Energy2 tb_;
  Energy2 ub_;

  // The diboson rapidity:
  double Yb_;
  // Note Yb_ = + lastY() if flipped = false 
  // but  Yb_ = - lastY() if flipped = true.
  // Yb_ is always defined with the quark travelling in the +z direction!

  // Masses of the final state bosons:
  Energy2 k12b_;
  Energy2 k22b_;

  // Polar and azimuthal angles of the dibosons in their rest frame:
  double theta1b_;
  double theta2b_;

public:

  /**
   * Default constructor 
   */
  born2to2Kinematics();

  /**
   * Meaningful constructor: takes the momenta 
   * and Bjorken x values of the partons involved
   * in the 2->2 scattering of the leading order
   * / virtual process and calculates all interesting
   * Mandelstam and Born variables.
   */
  born2to2Kinematics(vector<Lorentz5Momentum> Momenta, double xp, double xm);

public:
  /**
   * Read-only access to all of the above member variables.
   */

  // Leading order momentum fractions and associated etabar's:
  inline double xpb()   const { return xpb_; }
  inline double etapb() const { return etapb_; } 
  inline double xmb()   const { return xmb_; }
  inline double etamb() const { return etamb_; }

  // The Born momenta according to the notation of the FMNR papers,
  // in the diboson centre of mass frame:
  inline Lorentz5Momentum p1b() const { return p1b_; }
  inline Lorentz5Momentum p2b() const { return p2b_; }
  inline Lorentz5Momentum k1b() const { return k1b_; }
  inline Lorentz5Momentum k2b() const { return k2b_; }

  // The diboson invariant mass / shat, that and uhat:
  inline Energy2 sb() const { return sb_; }
  inline Energy2 tb() const { return tb_; }
  inline Energy2 ub() const { return ub_; }

  // The diboson rapidity:
  inline double Yb() const { return Yb_; }
  // Note Yb_ = + lastY() if flipped = false 
  // but  Yb_ = - lastY() if flipped = true.
  // Yb_ is always defined with the quark travelling in the +z direction!

  // Masses of the final state bosons:
  inline Energy2 k12b() const { return k12b_; }
  inline Energy2 k22b() const { return k22b_; }

  // Polar and azimuthal angles of the dibosons in their rest frame:
  inline double theta1b() const { return theta1b_; }
  inline double theta2b() const { return theta2b_; }

  // A check to make sure that the momenta calculated from 
  // the energies and angles are equal to those of meMomenta().
  void sanityCheck() const;
     
};



/** \ingroup NLO2to2Kinematics
 *  The Real2to3Kinematics class is used to store information on the 
 *  the kinematics of the real emission processes needed for the 
 *  evaluation of matrix elements in the real part of the NLO process.
 */
class real2to3Kinematics {
 
private:
  /**
   * The born2to2Kinematics object underlying the 2->3 kinematics.
   */
  born2to2Kinematics bornVariables_;

  /**
   * The lower bound on the x integration.
   */
  double xbar_;

  // The `raw' radiative variables.
  double xt_;
  double y_;
  double xr_;

  /**
   * The momentum fraction of the parton incident from the +z direction.
   */
  double xpr_;

  /**
   * The momentum fraction of the parton incident from the -z direction.
   */
  double xmr_;

  /**
   * Invariants required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */
  // First the Born variables:
  Energy2 s2r_;
  Energy2 k12r_;
  Energy2 k22r_;
  double  theta1r_;
  double  theta2r_;
  // Then the rest:
  Energy2 sr_;
  Energy2 tkr_;
  Energy2 ukr_;
  double  cpsir_;
  double  cpsiprr_;
  double  betaxr_;
  double  v1r_;
  double  v2r_;
  Energy2 q1r_;
  Energy2 q2r_;
  Energy2 q1hatr_;
  Energy2 q2hatr_;
  Energy2 w1r_;
  Energy2 w2r_;
  Lorentz5Momentum p1r_;
  Lorentz5Momentum p2r_;
  Lorentz5Momentum kr_;
  Lorentz5Momentum k1r_;
  Lorentz5Momentum k2r_;

public:

  /**
   * Default constructor 
   */
  real2to3Kinematics();

  /**
   * Meaningful constructor: takes the Born variables
   * from the leading order /virtual 2->2 process and 
   * the raw \tilde{x}, y radiative variables and turns
   * these into a set of 2->3 momenta with associated 
   * Mandelstam variables, Bjorken x values etc. 
   * @param xt The \tilde{x} radiative variable.
   * @param y  The angular radiative variable (the cosine 
   * of the  polar angle of the emitted gluon in the partonic 
   * CMS frame). 
   */
  real2to3Kinematics(born2to2Kinematics bornVariables,double xt, double y);

  // A check to make sure that the momenta calculated from 
  // the energies and angles are equal to those of meMomenta().
  void sanityCheck() const;
     
public:
  /**
   * Read-only access to all of the above member variables.
   */

  // The born2to2Kinematics underlying the 2->3 kinematics
  inline born2to2Kinematics bornVariables() const { return bornVariables_; }

  // The lower bound on the x integration:
  inline double xbar() const { return xbar_; }

  // The `raw' radiative variables.
  inline double xt() const { return xt_; }
  inline double y() const { return y_; }
  inline double xr() const { return xr_; }

  /**
   * The momentum fraction of the parton incident from the +z direction.
   */
  inline double xpr() const { return xpr_; }

  /**
   * The momentum fraction of the parton incident from the -z direction.
   */
  inline double xmr() const { return xmr_; }

  /**
   * Invariants required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */

  // First the Born variables:
  inline Energy2 s2r() const { return s2r_; }
  inline Energy2 k12r() const { return k12r_; }
  inline Energy2 k22r() const { return k22r_; }
  inline double  theta1r() const { return theta1r_; }
  inline double  theta2r() const { return theta2r_; }

  // Then the rest:
  inline Energy2 sr() const { return sr_; }
  inline Energy2 tkr() const { return tkr_; }   
  inline Energy2 ukr() const { return ukr_; }
  inline double  cpsir() const { return cpsir_; }
  inline double  cpsipr() const { return cpsiprr_; }
  inline double  betaxr() const { return betaxr_; }
  inline double  v1r() const { return v1r_; }
  inline double  v2r() const { return v2r_; }
  inline Energy2 q1r() const { return q1r_; }
  inline Energy2 q2r() const { return q2r_; } 
  inline Energy2 q1hatr() const { return q1hatr_; }
  inline Energy2 q2hatr() const { return q2hatr_; }
  inline Energy2 w1r() const { return w1r_; }
  inline Energy2 w2r() const { return w2r_; }
  inline Lorentz5Momentum p1r() const { return p1r_; }
  inline Lorentz5Momentum p2r() const { return p2r_; }
  inline Lorentz5Momentum kr()  const { return kr_ ; }
  inline Lorentz5Momentum k1r() const { return k1r_; }
  inline Lorentz5Momentum k2r() const { return k2r_; }
};

}

#endif /* HERWIG_NLO2to2Kinematics_H */
