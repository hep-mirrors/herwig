// -*- C++ -*-
#ifndef HERWIG_IS_QtildaShowerKinematics1to2_H
#define HERWIG_IS_QtildaShowerKinematics1to2_H
//
// This is the declaration of the <!id>IS_QtildaShowerKinematics1to2<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class provides the specific Intial State shower <BR>
// kinematics information.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:QtildaShowerKinematics1to2.html">QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:QtildaShowerKinematics1to2.html">FS_QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:KinematicsReconstructor.html">KinematicsReconstructor.h</a>.
// 

#include "QtildaShowerKinematics1to2.h"


namespace Herwig {

using namespace ThePEG;

class IS_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  inline IS_QtildaShowerKinematics1to2();
  inline IS_QtildaShowerKinematics1to2(const IS_QtildaShowerKinematics1to2 &);
  virtual ~IS_QtildaShowerKinematics1to2();
  // Standard ctors and dtor.

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

  virtual void updateParent( tCollecShoKinPtr & shoKinChildren );
  // update the parent Kinematics from the knowledge of the kinematics
  // of the children.  This method will be used by the 
  // <!class>KinematicsReconstructor<!!class>.


  virtual Energy jetMass();
  // The method returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> is true, 
  // and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

private:

  IS_QtildaShowerKinematics1to2 & operator=(const IS_QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

};

}

#include "IS_QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_IS_QtildaShowerKinematics1to2_H */
