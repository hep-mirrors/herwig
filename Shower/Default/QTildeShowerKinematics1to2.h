// -*- C++ -*-
//
// QTildeShowerKinematics1to2.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeShowerKinematics1to2_H
#define HERWIG_QTildeShowerKinematics1to2_H
//
// This is the declaration of the QTildeShowerKinematics1to2 class.
//

#include "Herwig/Shower/Base/ShowerKinematics.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "QTildeShowerKinematics1to2.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This abstract class describes the common features for initial and final
 *  state radiation kinematics for \f$1\to2\f$ branchings and for
 *  the choice of \f$\tilde{q}\f$ as evolution variable.
 *
 *  @see ShowerKinematics
 *  @see IS_QTildeShowerKinematics1to2
 *  @see FS_QTildeShowerKinematics1to2
 *  @see KinematicsReconstructor
 */ 
class QTildeShowerKinematics1to2: public ShowerKinematics {

public:

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
  const Lorentz5Momentum & pVector() const {return _pVector;}

  /**
   * Access to the \f$n\f$ vector used to describe the kinematics.
   */
  const Lorentz5Momentum & nVector() const {return _nVector;}

  /**
   *  Dot product of thew basis vectors
   */
  Energy2 p_dot_n() const {return _pVector*_nVector;}

  /**
   * Converts a Sudakov parametrization of a momentum w.r.t. the given 
   * basis \f$p\f$ and \f$n\f$ into a 5 momentum.
   * @param alpha The \f$\alpha\f$ parameter of the Sudakov parameterisation
   * @param beta  The \f$\beta\f$ parameter of the Sudakov parameterisation
   * @param px    The \f$x\f$-component of the transverse momentum in the Sudakov 
   *              parameterisation
   * @param py    The \f$x\f$-component of the transverse momentum in the Sudakov 
   *              parameterisation
   */
  Lorentz5Momentum sudakov2Momentum(double alpha, double beta,
				    Energy px, Energy py) const {
    return alpha*_pVector + beta*_nVector + px*_xPerp+py*_yPerp;
  }

  /**
   *  Transform the shower kinematics (usually the reference vectors)
   */
  virtual void transform(const LorentzRotation & r);

protected:

  /**
   *  Set the basis vectors
   */
  void setBasis(const Lorentz5Momentum &p, const Lorentz5Momentum & n,
		Frame frame);

  /**
   *  Set a preliminary momentum for the particle
   */
  void setMomentum(tShowerParticlePtr,bool timelike) const;

  /**
   *  Construct the spin info object for a shower particle
   */
  void constructSpinInfo(tShowerParticlePtr,bool timelike) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeShowerKinematics1to2 & operator=(const QTildeShowerKinematics1to2 &);

private:

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum _pVector;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum _nVector;

  /**
   *  x \f$q_\perp\f$ reference vector
   */
  LorentzVector<double> _xPerp;

  /**
   *  y \f$q_\perp\f$reference vector
   */
  LorentzVector<double> _yPerp;
};

}

#endif /* HERWIG_QTildeShowerKinematics1to2_H */
