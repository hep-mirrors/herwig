// -*- C++ -*-
//
// ShowerKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the ShowerKinematics class.
//

#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "ShowerKinematics.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This is the abstract base class from which all other shower
 * kinematics classes derive. The main purpose of the
 * shower kinematics classes is to allow the reconstruction
 * of jet masses, at the end of the showering (indeed, for
 * multi-scale showering, at the end of each scale-range evolution).
 * This is necessary for the kinematics reshuffling
 * in order to compensate the recoil of the emissions.
 * The KinematicsReconstructor class is in 
 * charge of this job, and which is the main "user" of
 * ShowerKinematics and its derived classes.
 * How this is done depends on the choice of kinematics variables 
 * and whether the jet is time-like (forward evolved) or 
 * space-like (backward evolved), whereas the class ShowerKinematics
 * describes only the common features which are independent by them.   
 *
 *  In general there are a number of methods specific to a shower approach
 *
 * @see KinematicsReconstructor
 */
class ShowerKinematics: public Base {

public:
  /**
   * The default constructor.
   */
  ShowerKinematics() : Base(), _isTheJetStartingPoint( false ),
		       _scale(), _z( 0.0 ), _phi( 0.0 ), _pt(),
		       _splitFun() {}

  /**
   *  The updateChildren and updateParent
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
			      const ShowerParticleVector & theChildren) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ShowerParticleVector & theChildren) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   * @param px The \f$x\f$ component of the \f$p_T\f$.
   * @param py The \f$y\f$ component of the \f$p_T\f$.
   */
  virtual void updateLast(const tShowerParticlePtr theLast,
			  Energy px, Energy py) const;
  //@}

  /**
   *  The reconstructLast, reconstructChildren and reconstructParent members
   *  are used during the reconstruction 
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
  virtual void reconstructChildren(const tShowerParticlePtr theParent, 
			      const ShowerParticleVector & theChildren) const;

  /**
   * Reconstruct the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void reconstructParent(const tShowerParticlePtr theParent, 
				 const ParticleVector & theChildren) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   * @param mass The mass to be used, if less than zero on-shell
   */
  virtual void reconstructLast(const tShowerParticlePtr theLast,
			       unsigned int iopt, Energy mass=-1.*GeV) const;

  /**
   *  Perform any initial calculations needed after the branching has been selected
   * @param particle The branching particle
   * @param parent The bema particle for the jet if needed
   */
  virtual void initialize(ShowerParticle & particle,PPtr parent);
  //@}

public:

  /**
   * Set/access the flag that tells whether or not this ShowerKinematics
   * object is associated to the starting particle of the jet: only in this
   * case it is sensible to use the two main virtual methods below.
   */
  //@{
  /**
   * Set the starting point flag
   */
  void isTheJetStartingPoint(const bool );
  
  /**
   * Get the starting point flag
   */
  bool isTheJetStartingPoint() const;
  //@}

  /**
   * Virtual function to return a set of basis vectors, specific to
   * the type of evolution. This function will be used by the
   * ForwardShowerEvolver in order to access \f$p\f$ and \f$n\f$, 
   * which in turn are members of the concrete class QtildaShowerKinematics1to2.
   */
  virtual vector<Lorentz5Momentum> getBasis() const = 0;


  /**
   *  Set/Get methods for the kinematic variables
   */
  //@{
  /**
   * Access the scale of the splitting.
   */
  Energy scale() const { return _scale; }

  /**
   * Set the scale of the splitting.
   */
  void scale(const Energy in) { _scale=in; }

  /**
   *  Access the energy fraction, \f$z\f$.
   */
  double z() const { return _z; }

  /**
   *  Set the energy fraction, \f$z\f$.
   */
  void z(const double in) { _z=in; }

  /**
   *  Access the azimuthal angle, \f$\phi\f$.
   */
  double phi() const { return _phi; }

  /**
   *  Set the azimuthal angle, \f$\phi\f$.
   */
  void phi(const double in) { _phi=in; }

  /**
   *  Access the relative \f$p_T\f$ for the branching
   */
  Energy pT() const { return _pt; }

  /**
   *  Set the relative \f$p_T\f$ for the branching
   */
  void pT(const Energy in) { _pt=in; }
  //@}

  /**
   *  Set and get methods for the SplittingFunction object
   */
  /**
   * Access the SplittingFunction object responsible of the 
   * eventual branching of this particle.
   */
  tSplittingFnPtr splittingFn() const { return _splitFun; }

  /**
   * Set the SplittingFunction object responsible of the 
   * eventual branching of this particle.
   */
  void splittingFn(const tSplittingFnPtr sf) { _splitFun=sf; }
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerKinematics & operator=(const ShowerKinematics &);

private:

  /**
   *  Is this the starting point of the jet
   */
  bool _isTheJetStartingPoint;

  /**
   *  The \f$\tilde{q}\f$ evolution variable.
   */
  Energy _scale;

  /**
   *  The energy fraction, \f$z\f$
   */
  double _z;

  /**
   *  The azimuthal angle, \f$\phi\f$.
   */
  double _phi;

  /**
   *  The relative \f$p_T\f$
   */
  Energy _pt;

  /**
   *  The splitting function for the branching of the particle
   */
  tSplittingFnPtr _splitFun;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerKinematics. */
template <>
struct BaseClassTrait<Herwig::ShowerKinematics,1> {
  /** Typedef of the first base class of ShowerKinematics. */
  typedef Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerKinematics class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerKinematics>
  : public ClassTraitsBase<Herwig::ShowerKinematics> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerKinematics"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerKinematics is implemented. It may also include several, space-separated,
   * libraries if the class ShowerKinematics depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ShowerKinematics_H */
