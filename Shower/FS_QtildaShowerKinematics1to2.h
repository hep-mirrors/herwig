// -*- C++ -*-
#ifndef HERWIG_FS_QtildaShowerKinematics1to2_H
#define HERWIG_FS_QtildaShowerKinematics1to2_H
//
// This is the declaration of the <!id>FS_QtildaShowerKinematics1to2<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the specific Final State shower <BR>
// kinematics information.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:QtildaShowerKinematics1to2.html">QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:QtildaShowerKinematics1to2.html">IS_QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:KinematicsReconstructor.html">KinematicsReconstructor.h</a>.
// 

#include "QtildaShowerKinematics1to2.h"


namespace Herwig {

using namespace Pythia7;

class FS_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  inline FS_QtildaShowerKinematics1to2();
  inline FS_QtildaShowerKinematics1to2(const FS_QtildaShowerKinematics1to2 &);
  virtual ~FS_QtildaShowerKinematics1to2();
  // Standard ctors and dtor.

  inline FS_QtildaShowerKinematics1to2( const Lorentz5Momentum & p, 
					const Lorentz5Momentum & n );
  // Creator with the two defining vectors <!id>p<!!id> and <!id>n<!!id> 

  virtual void updateChildren( const double parentSudAlpha, 
			       const Energy parentSudPx, const Energy parentSudPy, 
                               vector<double> & sudAlphaVect, 
			       vector<Energy> & sudPxVect, vector<Energy> & sudPyVect );
  // Along with the showering evolution --- going forward for
  // time-like (forward) evolution, and going backward for space-like
  // (backward) evolution --- the Sudakov variables associated to the
  // branching products are calcalted and returned, from the knowledge
  // of the parent Sudakov variables.   
  // Note that only <I>alpha</I> and <I>p_perp</I> can be reconstructed 
  // at this moment and we will obtain instead beta only later, 
  // using <!id>updateParent()<!!id>.

  virtual void updateChildren( const tShowerParticlePtr theParent, 
			       const ParticleVector theChildren );
  // Along with the showering evolution --- going forward for
  // time-like (forward) evolution, and going backward for space-like
  // (backward) evolution --- the kinematical variables of the
  // branching products are calculated and updated from the knowledge
  // of the parent kinematics.  This method is used by the
  // <!class>ForwardShowerEvolver<!!class>.  ***ACHTUNG*** Might be
  // extended to update colour connections as well.

  virtual void updateParent( const tShowerParticlePtr theParent, 
			     const ParticleVector theChildren );
  // update the parent Kinematics from the knowledge of the kinematics
  // of the children.  This method will be used by the 
  // <!class>KinematicsReconstructor<!!class>.

  virtual void updateLast( const tShowerParticlePtr theLast );
  // update the kinematical data of a particle when a reconstruction
  // fixpoint was found.  This will highly depend on the kind of
  // kinematics chosen and will be defined in the inherited concrete
  // classes. This method will be used by the
  // <!class>KinematicsReconstructor<!!class>.

  virtual Energy jetMass();
  // The method returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> 
  // is true, and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

private:

  FS_QtildaShowerKinematics1to2 & operator=(const FS_QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

};

}

#include "FS_QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_FS_QtildaShowerKinematics1to2_H */
