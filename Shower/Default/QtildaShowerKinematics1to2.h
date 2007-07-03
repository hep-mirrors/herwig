// -*- C++ -*-
#ifndef HERWIG_QtildaShowerKinematics1to2_H
#define HERWIG_QtildaShowerKinematics1to2_H
//
// This is the declaration of the QtildaShowerKinematics1to2 class.
//

#include "Herwig++/Shower/Base/ShowerKinematics.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "QtildaShowerKinematics1to2.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This abstract class describes the common features for initial and final
 *  state radiation kinematics for \f$1\to2\f$ branchings and for
 *  the choice of \f$\tilde{q}\f$ as evolution variable.
 *
 *  @see ShowerKinematics
 *  @see IS_QtildaShowerKinematics1to2
 *  @see FS_QtildaShowerKinematics1to2
 *  @see KinematicsReconstructor
 */ 
class QtildaShowerKinematics1to2: public ShowerKinematics {

public:

  /**
   * Default Constructor
   */
  inline QtildaShowerKinematics1to2();

  /**
   * Implementation of the virtual function returning a set of basis vectors, specific to
   * the type of evolution.  This function will be used by the
   * ForwardShowerEvolver in order to access \f$p\f$
   * and \f$n\f$.
   */
  virtual vector<Lorentz5Momentum> getBasis() const; 

  /**
   * Access to the \f$p\f$ vector used to describe the kinematics.
   */
  inline const Lorentz5Momentum & pVector() const;

  /**
   * Access to the \f$n\f$ vector used to describe the kinematics.
   */
  inline const Lorentz5Momentum & nVector() const;

  /**
   *  Dot product of thew basis vectors
   */
  inline const Energy2 p_dot_n() const;

  /**
   * Converts a Sudakov parametrization of a momentum w.r.t. the given 
   * basis \f$p\f$ and \f$n\f$ into a 5 momentum.
   * @param alpha The \f$\alpha\f$ parameter of the Sudakov parameterisation
   * @param beta  The \f$\beta\f$ parameter of the Sudakov parameterisation
   * @param px    The \f$x\f$-component of the transverse momentum in the Sudakov 
   *              parameterisation
   * @param py    The \f$x\f$-component of the transverse momentum in the Sudakov 
   *              parameterisation
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   */
  Lorentz5Momentum sudakov2Momentum(double alpha, double beta, Energy px, Energy py,
				    unsigned int iopt) const;

protected:

  /**
   *  Set the basis vectors
   */
  void setBasis(const Lorentz5Momentum &p, const Lorentz5Momentum & n);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QtildaShowerKinematics1to2 & operator=(const QtildaShowerKinematics1to2 &);

private:

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum _pVector;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum _nVector;

};

}

#include "QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
