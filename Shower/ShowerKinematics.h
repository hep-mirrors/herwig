// -*- C++ -*-
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the ShowerKinematics class.

#include "ShowerConfig.h"
#include "ThePEG/Pointer/Ptr.h"
#include "ThePEG/Pointer/ReferenceCounted.h"
#include "ThePEG/Pointer/PtrTraits.h"
#include "ThePEG/Pointer/RCPtr.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ShowerParticle.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * This is the abstract base class from which all other shower
 * kinematics classes derive from. The main purpose of the
 * shower kinematics classes is to allow the reconstruction
 * of jet masses, at the end of the showering (indeed, for
 * multi-scale showering, at the end of each scale-range evolution). <BR> 
 * --- This is necessary for the kinematics reshuffling
 * in order to compensate the recoil of the emissions.
 * It is the class KinematicsReconstructor which is in 
 * charge of this job, and which is the main "user" of
 * ShowerKinematics and its derived classes. --- <BR>
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

  /**
   * Standard ctors and dtor.
   */
  inline ShowerKinematics();
  inline ShowerKinematics(const ShowerKinematics &);
  virtual ~ShowerKinematics();

  /**
   * Creator with the two defining vectors p and n. 
   */
  inline ShowerKinematics(const Lorentz5Momentum & p, 
			  const Lorentz5Momentum & n);

  /**
   * Set/access the flag that tells whether or not this ShowerKinematics
   * object is associated to the starting particle of the jet: only in this
   * case it is sensible to use the two main virtual methods below.
   */
  void isTheJetStartingPoint(const bool inputIsTheJetStartingPoint);
  bool isTheJetStartingPoint() const;

  virtual void updateChildren(const double parentAlpha, const Energy parentPx, 
			      const Energy parentPy, 
			      vector<double> & alphaVect, 
			      vector<Energy> & pxVect, 
			      vector<Energy> & pyVect) = 0;

  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the Sudakov variables associated to the
   * branching products are calcalted and returned, from the knowledge
   * of the parent Sudakov variables.   
   * Note that only <I>alpha</I> and <I>p_perp</I> can be reconstructed 
   * at this moment and we will obtain instead beta only later, 
   * using updateParent().
   */

  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics.  This method is used by the
   * ForwardShowerEvolver.  ***ACHTUNG*** Might be
   * extended to update colour connections as well.
   */
  virtual void updateChildren(const tShowerParticlePtr theParent, 
			      const ParticleVector theChildren) = 0;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   */
  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ParticleVector theChildren) = 0;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
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
   * ForwardShowerEvolver in order to access <I>p</I> and <I>n</I>, 
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
   */
  Lorentz5Momentum referenceFrame(const Lorentz5Momentum & particleMomentum,
				  const vector<Lorentz5Momentum> & partnerP,
				  const vector<Energy> & evolutionScales);

  /**
   * Access to/set the scale of the splitting.
   */
  inline Energy qtilde() const;
  inline void qtilde( const Energy inputQtilde );

  /**
   * Access to/set the resolution scale of the splitting.
   */
  inline Energy resScale() const;
  inline void setResScale( const Energy inputQ );

  /**
   * Access to/set a cutoff scale of the splitting.
   */
  inline Energy kinScale() const;
  inline void setKinScale( const Energy inputQ );

  inline double z() const;
  inline void z( const double inputZ );
  inline double phi() const;
  inline void phi( const double inputPhi );

private:

  /**
   * Private and non-existent assignment operator.
   */
  ShowerKinematics & operator=(const ShowerKinematics &);

  bool _isTheJetStartingPoint;

  Energy _qtilde;
  Energy _q0;
  Energy _kinQ0;
  double _z;
  double _phi;
};

}

#include "ShowerKinematics.icc"

#endif /* HERWIG_ShowerKinematics_H */
