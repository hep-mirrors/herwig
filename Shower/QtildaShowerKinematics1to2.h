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
#include "Pythia7/CLHEPWrap/Lorentz5Vector.h"


namespace Herwig {

using namespace Pythia7;

class QtildaShowerKinematics1to2: public ShowerKinematics {

public:

  inline QtildaShowerKinematics1to2();
  inline QtildaShowerKinematics1to2(const QtildaShowerKinematics1to2 &);
  virtual ~QtildaShowerKinematics1to2();
  // Standard ctors and dtor.

  inline QtildaShowerKinematics1to2( const Lorentz5Momentum & p, 
				     const Lorentz5Momentum & n,
				     const Energy inputOnShellMass );
  // Creator with the two defining vectors <!id>p<!!id> and <!id>n<!!id> 
  // and the on-shell mass of the particle.

  virtual void updateChildren( const double parentSudAlpha, 
			       const Energy parentSudPx, const Energy parentSudPy, 
                               vector<double> & sudAlphaVect, 
			       vector<Energy> & sudPxVect, vector<Energy> & sudPyVect ) = 0;
  // Along with the showering evolution --- going forward for
  // time-like (forward) evolution, and going backward for space-like
  // (backward) evolution --- the Sudakov variables associated to the
  // branching products are calcalted and returned, from the knowledge
  // of the parent Sudakov variables.   
  // Note that only <I>alpha</I> and <I>p_perp</I> can be reconstructed 
  // at this moment and we will obtain instead beta only later, 
  // using <!id>updateParent()<!!id>.

  virtual void updateParent( tCollecShoKinPtr & shoKinChildren ) = 0;
  // update the parent Kinematics from the knowledge of the kinematics
  // of the children.  This method will be used by the 
  // <!class>KinematicsReconstructor<!!class>.

  virtual Energy jetMass() = 0;
  // Pure virtual method, to be defined in a derived class.
  // The method returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> is true, 
  // and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

  inline double z() const;
  inline void z( const double inputZ );
  inline double phi() const;
  inline void phi( const double inputPhi );
  // Access/set to the generated kinematics variables of the splitting <I>1-&GT;2</I>.

  inline const Lorentz5Momentum & pVector() const;
  inline const Lorentz5Momentum & nVector() const;
  // Access to the <!id>p<!!id> and <!id>n<!!id> vectors used to described the kinematics.

private:

  QtildaShowerKinematics1to2 & operator=(const QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

  double _z;
  double _phi;
  Lorentz5Momentum _pVector;  //***LOOKHERE*** Re-think where to put them. 
  Lorentz5Momentum _nVector;  

};

}

#include "QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
