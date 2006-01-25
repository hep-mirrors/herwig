// -*- C++ -*-
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the ShowerKinematics class.
//

#include "ShowerConfig.h"
#include "ThePEG/Pointer/Ptr.h"
#include "ThePEG/Pointer/ReferenceCounted.h"
#include "ThePEG/Pointer/PtrTraits.h"
#include "ThePEG/Pointer/RCPtr.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ShowerParticle.h"
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
 * @see QtildaShowerKinematics1to2
 * @see KinematicsReconstructor
 */
class ShowerKinematics: public ReferenceCounted {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerKinematics();

  /**
   * The copy constructor.
   */
  inline ShowerKinematics(const ShowerKinematics &);

  /**
   * The destructor.
   */
  virtual ~ShowerKinematics();
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
   *  The updateChildren members to update the values of the \f$\alpha\f$ and 
   *  \f$p_\perp\f$ variables during the shower evolution.
   */
  //@{
  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the Sudakov variables associated to the
   * branching products are calculated and returned, from the knowledge
   * of the parent Sudakov variables.   
   * Note that only \f$\alpha\f$ and \f$p_{\perp}\f$ can be reconstructed 
   * at this moment and we will obtain \f$\beta\f$ later, 
   * using updateParent().
   * @param parentAlpha  The \f$\alpha\f$ of the parent
   * @param parentPx     The x-component of the \f$p_\perp\f$ of the parent
   * @param parentPy     The y-component of the \f$p_\perp\f$ of the parent
   * @param alphaVect    The \f$\alpha\f$'s of the children
   * @param pxVect       The x-components of the \f$p_\perp\f$ of the children 
   * @param pyVect       The y-components of the \f$p_\perp\f$ of the children 
   */
  virtual void updateChildren(const double parentAlpha, const Energy parentPx, 
			      const Energy parentPy, 
			      vector<double> & alphaVect, 
			      vector<Energy> & pxVect, 
			      vector<Energy> & pyVect) = 0;

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
  //@}

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

  /**
   * Pure virtual method to be defined in derived classes.
   * It returns the mass of jet. 
   * It should be used only if isTheJetStartingPoint() is true, 
   * and only after the jet kinematics reconstruction.
   * (performed by the KinematicsReconstructor class).
   */
  virtual Energy jetMass() = 0;

  /**
   * Virtual function to return a set of basis vectors, specific to
   * the type of evolution. This function will be used by the
   * ForwardShowerEvolver in order to access \f$p\f$ and \f$n\f$, 
   * which in turn are members of the concrete class QtildaShowerKinematics1to2.
   */
  virtual vector<Lorentz5Momentum> getBasis() = 0;

  /**
   * It returns the 5-vector momentum of the CM reference frame,
   * with respect to the Lab frame, in which the jet evolution
   * is supposed to be described. In order to determine such frame
   * it needs to know the momentum of the particle (to which this
   * ShowerKinematics object is attached to; the particle should be
   * the originator of the jet), and the momenta of the partners and 
   * the corresponding evolution scales for each type of interaction
   * considered (QCD, QED, EWK, ...) as defined in ShowerIndex. 
   * This method should be used only if isTheJetStartingPoint() is true.
   * @param particleMomentum The momentum of the particle
   * @param partnerP         The momentum of the partner
   * @param evolutionScales  The evolution scales
   */
  Lorentz5Momentum referenceFrame(const Lorentz5Momentum & particleMomentum,
				  const vector<Lorentz5Momentum> & partnerP,
				  const vector<Energy> & evolutionScales);

  /**
   *  Set/Get methods for the kinematic variables
   */
  //@{
  /**
   * Access the scale of the splitting.
   */
  inline Energy qtilde() const;

  /**
   * Set the scale of the splitting.
   */
  inline void qtilde(const Energy);

  /**
   * Access the resolution scale of the splitting.
   */
  inline Energy resScale() const;

  /**
   * Set the resolution scale of the splitting.
   */
  inline void setResScale( const Energy inputQ );

  /**
   * Access the cutoff scale of the splitting.
   */
  inline Energy kinScale() const;

  /**
   * Set the cutoff scale of the splitting.
   */
  inline void setKinScale( const Energy inputQ );

  /**
   *  Access the energy fraction, \f$z\f$.
   */
  inline double z() const;

  /**
   *  Set the energy fraction, \f$z\f$.
   */
  inline void z(const double);
  /**
   *  Access the azimuthal angle, \f$\phi\f$.
   */
  inline double phi() const;

  /**
   *  Set the azimuthal angle, \f$\phi\f$.
   */
  inline void phi(const double);
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
  Energy _qtilde;

  /**
   *  The resolution scale
   */
  Energy _q0;

  /**
   *  The kinematic cut-off
   */
  Energy _kinQ0;

  /**
   *  The energy fraction, \f$z\f$
   */
  double _z;

  /**
   *  The azimuthal angle, \f$\phi\f$.
   */
  double _phi;
};

}

#include "ShowerKinematics.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerKinematics.tcc"
#endif

#endif /* HERWIG_ShowerKinematics_H */
