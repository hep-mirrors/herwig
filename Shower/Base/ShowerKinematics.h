// -*- C++ -*-
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the ShowerKinematics class.
//

#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Shower/ShowerVariables.h"
#include "ThePEG/Config/ThePEG.h"
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
 * @see KinematicsReconstructor
 */
class ShowerKinematics: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerKinematics(ShowerVarsPtr);
  //@}

public:

  /**
   *  The updateChildren, updateLast and updateParent
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
			      const ShowerParticleVector theChildren) const = 0;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ParticleVector theChildren) const = 0;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   */
  virtual void updateLast(const tShowerParticlePtr theLast, unsigned int iopt) const = 0;
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
   *  Access to the ShowerVariables
   */
  inline ShowerVarsPtr showerVariables() const;

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

  /**
   *  Access the relative \f$p_T\f$ for the branching
   */
  inline Energy pT() const;

  /**
   *  Set the relative \f$p_T\f$ for the branching
   */
  inline void pT(const Energy);
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

  /**
   *  The relative \f$p_T\f$
   */
  Energy _pt;

  /**
   *  Pointer to the ShowerVariables to allow access
   */
  ShowerVarsPtr _variables;
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
  typedef ReferenceCounted NthBase;
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
  static string library() { return "ShowerKinematics.so"; }
};

/** @endcond */

}

#include "ShowerKinematics.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerKinematics.tcc"
#endif

#endif /* HERWIG_ShowerKinematics_H */
