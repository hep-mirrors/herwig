// -*- C++ -*-
#ifndef HERWIG_QtildaShowerKinematics1to2_H
#define HERWIG_QtildaShowerKinematics1to2_H
//
// This is the declaration of the <!id>QtildaShowerKinematics1to2<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class describes the common features for initial and final <BR>
// state radiation kinematics for <I>1 -&GT; 2</I> branchings and for <BR>
// the choice of Qtilda as evolution variable. <BR>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerKinematics.html">ShowerKinematics.h</a>, <BR>
// <a href="http:IS_QtildaShowerKinematics1to2.html">IS_QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:FS_QtildaShowerKinematics1to2.html">FS_QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:KinematicsReconstructor.html">KinematicsReconstructor.h</a>.
// 

#include "ShowerKinematics.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {

using namespace ThePEG;

class QtildaShowerKinematics1to2: public ShowerKinematics {

public:

  inline QtildaShowerKinematics1to2();
  inline QtildaShowerKinematics1to2(const QtildaShowerKinematics1to2 &);
  virtual ~QtildaShowerKinematics1to2();
  // Standard ctors and dtor.

  inline QtildaShowerKinematics1to2(const Lorentz5Momentum & p, 
				    const Lorentz5Momentum & n);
  // Creator with the two defining vectors <!id>p<!!id> and <!id>n<!!id> 

  virtual void updateChildren(const double parentAlpha, 
			      const Energy parentPx, const Energy parentPy, 
                              vector<double> & alphaVect, 
			      vector<Energy> & pxVect, 
			      vector<Energy> & pyVect) = 0;
  // Along with the showering evolution --- going forward for
  // time-like (forward) evolution, and going backward for space-like
  // (backward) evolution --- the Sudakov variables associated to the
  // branching products are calcalted and returned, from the knowledge
  // of the parent Sudakov variables.   
  // Note that only <I>alpha</I> and <I>p_perp</I> can be reconstructed 
  // at this moment and we will obtain instead beta only later, 
  // using <!id>updateParent()<!!id>.

  virtual void updateChildren(const tShowerParticlePtr theParent, 
			      const ParticleVector theChildren) = 0;
  // Along with the showering evolution --- going forward for
  // time-like (forward) evolution, and going backward for space-like
  // (backward) evolution --- the kinematical variables of the
  // branching products are calculated and updated from the knowledge
  // of the parent kinematics.  This method is used by the
  // <!class>ForwardShowerEvolver<!!class>.  ***ACHTUNG*** Might be
  // extended to update colour connections as well.

  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ParticleVector theChildren) = 0;
  // update the parent Kinematics from the knowledge of the kinematics
  // of the children.  This method will be used by the 
  // <!class>KinematicsReconstructor<!!class>.

  virtual void updateLast(const tShowerParticlePtr theLast) = 0;
  // update the kinematical data of a particle when a reconstruction
  // fixpoint was found.  This will highly depend on the kind of
  // kinematics chosen and will be defined in the inherited concrete
  // classes. This method will be used by the
  // <!class>KinematicsReconstructor<!!class>.

  virtual Energy jetMass() = 0;
  // Pure virtual method, to be defined in a derived class.
  // The method returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> is true, 
  // and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

  virtual vector<Lorentz5Momentum> getBasis(); 
  // virtual function to return a set of basis vectors, specific to
  // the type of evolution.  this function will be used by the
  // <!class>ForwardShowerEvolver<!!class> in order to access <I>p</I>
  // and <I>n</I>, which in turn are members of the concrete class
  // <!class>QtildaShowerKinematics1to2<!!class>.

  // Access/set to the generated kinematics variables of the splitting <I>1-&GT;2</I>.

  inline const Lorentz5Momentum & pVector() const;
  inline const Lorentz5Momentum & nVector() const;
  inline const double p_dot_n() const;
  // Access to the <!id>p<!!id> and <!id>n<!!id> vectors used to describe the kinematics.

private:

  QtildaShowerKinematics1to2 & operator=(const QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

protected: 

  double _z;
  double _phi;
  const Lorentz5Momentum _pVector;
  const Lorentz5Momentum _nVector;

  Lorentz5Momentum sudakov2Momentum(double alpha, double beta, Energy px, Energy py);
  // converts a Sudakov parametrization of a momentum w.r.t. the given 
  // basis <!id>p<!!id> and <!id>n<!!id> into a 5 momentum.

};

}

#include "QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
