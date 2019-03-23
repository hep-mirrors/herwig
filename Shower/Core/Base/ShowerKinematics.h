// -*- C++ -*-
//
// ShowerKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the ShowerKinematics class.
//

#include "Herwig/Shower/Core/ShowerConfig.h"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig/Shower/Core/Base/SudakovFormFactor.h"
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
		       _sudakov() {}

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
   * @param parent   The parent
   * @param children The children
   * @param partnerType The type of evolution partner
   */
  virtual void updateChildren(const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children,
			      ShowerPartnerType partnerType,
			      bool massVeto ) const;

  virtual void resetChildren( const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param parent   The parent
   * @param children The children
   * @param partnerType The type of evolution partner
   */
  virtual void updateParent(const tShowerParticlePtr parent,
			    const ShowerParticleVector & children,
			    ShowerPartnerType partnerType) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param last The particle.
   * @param px The \f$x\f$ component of the \f$p_T\f$.
   * @param py The \f$y\f$ component of the \f$p_T\f$.
   */
  virtual void updateLast(const tShowerParticlePtr last,
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
   * @param parent   The parent
   * @param children The children
   */
  virtual void reconstructChildren(const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children) const;

  /**
   * Reconstruct the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param parent   The parent
   * @param children The children
   */
  virtual void reconstructParent(const tShowerParticlePtr parent, 
				 const ParticleVector & children) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param last The particle.
   * @param mass The mass to be used, if less than zero on-shell
   */
  virtual void reconstructLast(const tShowerParticlePtr last, Energy mass=-1.*GeV) const;
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
  void pT(const Energy in) const { _pt=in; }
  //@}

  /**
   *  Set and get methods for the SplittingFunction object
   */
  //@{
  /**
   * Access the SplittingFunction object responsible of the 
   * eventual branching of this particle.
   */
  tSplittingFnPtr splittingFn() const { return _sudakov-> splittingFn(); }
  //@}

  /**
   *  Set and get methods for the SudakovFormFactor object
   */
  /**
   * Access the SudakovFormFactor object responsible of the 
   * eventual branching of this particle.
   */
  tSudakovPtr SudakovFormFactor() const { return _sudakov; }

  /**
   * Set the SudakovFormFactor object responsible of the 
   * eventual branching of this particle.
   */
  void SudakovFormFactor(const tSudakovPtr sud) { _sudakov=sud; }
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerKinematics & operator=(const ShowerKinematics &) = delete;

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
  mutable Energy _pt;

  /**
   *  The splitting function for the branching of the particle
   */
  tSudakovPtr _sudakov;

};

}

#endif /* HERWIG_ShowerKinematics_H */
