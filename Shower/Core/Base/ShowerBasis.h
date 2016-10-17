// -*- C++ -*-
#ifndef Herwig_ShowerBasis_H
#define Herwig_ShowerBasis_H
//
// This is the declaration of the ShowerBasis class.
//

#include "ShowerBasis.fh"
#include "ShowerParticle.fh"
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ShowerBasis class stores the basis vectors used by the shower
 */
class ShowerBasis: public Base {

public:

  /**
   *   enum for the frame definition
   */
  enum Frame {BackToBack,Rest};

public:

  /**
   * The default constructor.
   */
  ShowerBasis() {}
  
  /**
   * Access to the frame definition
   */
  Frame frame() const {return frame_;}

  /**
   * Implementation of the virtual function returning a set of basis vectors, specific to
   * the type of evolution.  This function will be used by the
   * ForwardShowerEvolver in order to access \f$p\f$
   * and \f$n\f$.
   */
  virtual vector<Lorentz5Momentum> getBasis() const; 

  /**
   *  Set the basis vectors
   */
  void setBasis(const Lorentz5Momentum &p, const Lorentz5Momentum & n,
		Frame frame);

  /**
   * Access to the \f$p\f$ vector used to describe the kinematics.
   */
  const Lorentz5Momentum & pVector() const {return pVector_;}

  /**
   * Access to the \f$n\f$ vector used to describe the kinematics.
   */
  const Lorentz5Momentum & nVector() const {return nVector_;}

  /**
   *  Dot product of thew basis vectors
   */
  Energy2 p_dot_n() const {return pVector_*nVector_;}

  /**
   *  Transform the shower kinematics (usually the reference vectors)
   */
  virtual void transform(const LorentzRotation & r);

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
    return alpha*pVector_ + beta*nVector_ + px*xPerp_ + py*yPerp_;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerBasis & operator=(const ShowerBasis &);

private:

  /**
   *  The frame in which the basis vectors are defined
   */
  Frame frame_;

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum pVector_;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum nVector_;

  /**
   *  x \f$q_\perp\f$ reference vector
   */
  LorentzVector<double> xPerp_;

  /**
   *  y \f$q_\perp\f$reference vector
   */
  LorentzVector<double> yPerp_;

};

}

#endif /* Herwig_ShowerBasis_H */
