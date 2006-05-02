// -*- C++ -*-
#ifndef HERWIG_QtildaShowerKinematics1to2_H
#define HERWIG_QtildaShowerKinematics1to2_H
//
// This is the declaration of the QtildaShowerKinematics1to2 class.
//

#include "ShowerKinematics.h"
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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline QtildaShowerKinematics1to2(ShowerVarsPtr);

  /**
   * The copy constructor.
   */
  inline QtildaShowerKinematics1to2(const QtildaShowerKinematics1to2 &);

  /**
   * Constructor with the two defining vectors \f$p\f$ and \f$n\f$.
   * @param p The \f$p\f$ reference vector.
   * @param n The \f$n\f$ reference vector.
   * @param var Pointer to the ShowerVariables
   */
  inline QtildaShowerKinematics1to2(const Lorentz5Momentum & p, 
				    const Lorentz5Momentum & n,
				    ShowerVarsPtr var);

  /**
   * The destructor.
   */
  virtual ~QtildaShowerKinematics1to2();
  //@}

public:

  /**
   *  The updateChildren, updateParent and updateLast
   *  members to update the values of the \f$\alpha\f$ and 
   *  \f$p_\perp\f$ variables during the shower evolution.
   */
  //@{
  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics. 
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateChildren(const tShowerParticlePtr theParent, 
			      const ShowerParticleVector theChildren) = 0;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ParticleVector theChildren) = 0;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   */
  virtual void updateLast(const tShowerParticlePtr theLast) = 0;
  //@}

  /**
   * Implementation of the virtual function returning a set of basis vectors, specific to
   * the type of evolution.  This function will be used by the
   * ForwardShowerEvolver in order to access \f$p\f$
   * and \f$n\f$.
   */
  virtual vector<Lorentz5Momentum> getBasis(); 

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
  inline const double p_dot_n() const;

protected:

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
  Lorentz5Momentum sudakov2Momentum(double alpha, double beta, Energy px, Energy py);

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
  const Lorentz5Momentum _pVector;

  /**
   *  The \f$n\f$ reference vector
   */
  const Lorentz5Momentum _nVector;

};

}

#include "QtildaShowerKinematics1to2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtildaShowerKinematics1to2.tcc"
#endif

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
