// -*- C++ -*-
#ifndef HERWIG_QtildaShowerKinematics1to2_H
#define HERWIG_QtildaShowerKinematics1to2_H
//
// This is the declaration of the QtildaShowerKinematics1to2 class.

#include "ShowerKinematics.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This abstract class describes the common features for initial and final
 *  state radiation kinematics for <I>1 -&GT; 2</I> branchings and for
 *  the choice of Qtilda as evolution variable.
 *
 *  @see ShowerKinematics
 *  @see IS_QtildaShowerKinematics1to2
 *  @see FS_QtildaShowerKinematics1to2
 *  @see KinematicsReconstructor
 */ 
class QtildaShowerKinematics1to2: public ShowerKinematics {

public:

  /**
   * Standard ctors and dtor.
   */
  inline QtildaShowerKinematics1to2();
  inline QtildaShowerKinematics1to2(const QtildaShowerKinematics1to2 &);
  virtual ~QtildaShowerKinematics1to2();

  /**
   * Creator with the two defining vectors p and n.
   */
  inline QtildaShowerKinematics1to2(const Lorentz5Momentum & p, 
				    const Lorentz5Momentum & n);

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
  virtual void updateChildren(const double parentAlpha, 
			      const Energy parentPx, const Energy parentPy, 
                              vector<double> & alphaVect, 
			      vector<Energy> & pxVect, 
			      vector<Energy> & pyVect) = 0;

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
   * fixpoint was found.  This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by theKinematicsReconstructor.
   */
  virtual void updateLast(const tShowerParticlePtr theLast) = 0;

  /**
   * Pure virtual method, to be defined in a derived class.
   * The method returns the mass of jet. 
   * It should be used only if isTheJetStartingPoint() is true, 
   * and only after the jet kinematics reconstruction.
   * (performed by the KinematicsReconstructor class).
   */
  virtual Energy jetMass() = 0;

  /**
   * Virtual function to return a set of basis vectors, specific to
   * the type of evolution.  this function will be used by the
   * ForwardShowerEvolver in order to access <I>p</I>
   * and <I>n</I>, which in turn are members of the concrete class
   * QtildaShowerKinematics1to2.
   */
  virtual vector<Lorentz5Momentum> getBasis(); 

  /**
   * Access/set to the generated kinematics variables of the splitting <I>1-&GT;2</I>.
   */

  /**
   * Access to the p and n vectors used to describe the kinematics.
   */
  inline const Lorentz5Momentum & pVector() const;
  inline const Lorentz5Momentum & nVector() const;
  inline const double p_dot_n() const;

private:

  /**
   * Private and non-existent assignment operator.
   */
  QtildaShowerKinematics1to2 & operator=(const QtildaShowerKinematics1to2 &);

protected: 

  double _z;
  double _phi;
  const Lorentz5Momentum _pVector;
  const Lorentz5Momentum _nVector;

  /**
   * Converts a Sudakov parametrization of a momentum w.r.t. the given 
   * basis p and n into a 5 momentum.
   */
  Lorentz5Momentum sudakov2Momentum(double alpha, double beta, Energy px, Energy py);

};

}

#include "QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
