// -*- C++ -*-
#ifndef HERWIG_ShowerKinematics_H
#define HERWIG_ShowerKinematics_H
//
// This is the declaration of the <!id>ShowerKinematics<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is the abstract base class from which all other shower <BR>
// kinematics classes derive from. The main purpose of the <BR>
// shower kinematics classes is to allow the reconstruction <BR>
// of jet masses, at the end of the showering (indeed, for <BR>
// multi-scale showering, at the end of each scale-range evolution). <BR> 
// --- This is necessary for the kinematics reshuffling <BR>
// in order to compensate the recoil of the emissions. <BR>
// It is the class <!class>KinematicsReconstructor<!!class> which is in <BR> 
// charge of this job, and which is the main "user" of <BR>
// <!id>ShowerKinematics<!!id> and its derived classes. --- <BR>
// How this is done depends on the choice of kinematics variables <BR>
// and whether the jet is time-like (forward evolved) or <BR>
// space-like (backward evolved), whereas the class <!id>ShowerKinematics<!!id> <BR>
// describes only the common features which are independent by them.   
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:QtildaShowerKinematics1to2.html">QtildaShowerKinematics1to2.h</a>, <BR>
// <a href="http:KinematicsReconstructor.html">KinematicsReconstructor.h</a>.
// 

#include "ShowerConfig.h"
#include "Pythia7/Pointer/Ptr.h"
#include "Pythia7/Pointer/ReferenceCounted.h"
#include "Pythia7/Pointer/PtrTraits.h"
#include "Pythia7/Pointer/RCPtr.h"
#include "Herwig++/Config/GlobalParameters.h"


namespace Herwig {

using namespace Pythia7;

class ShowerKinematics: public ReferenceCounted {

public:

  inline ShowerKinematics();
  inline ShowerKinematics(const ShowerKinematics &);
  virtual ~ShowerKinematics();
  // Standard ctors and dtor.

  void isTheJetStartingPoint(const bool inputIsTheJetStartingPoint);
  bool isTheJetStartingPoint() const;
  // Set/access the flag that tells whether or not this <!id>ShowerKinematics<!!id>
  // object is associated to the starting particle of the jet: only in this
  // case it is sensible to use the two main virtual methods below.

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
  // Pure virtual method to be defined in derived classes.
  // It returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> is true, 
  // and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

  Lorentz5Momentum referenceFrame( const Lorentz5Momentum & particleMomentum,
				   const vector<Lorentz5Momentum> & partnersMomenta,
				   const vector<Energy> & evolutionScales );
  // It returns the 5-vector momentum of the CM reference frame,
  // with respect to the Lab frame, in which the jet evolution
  // is supposed to be described. In order to determine such frame
  // it needs to know the momentum of the particle (to which this
  // <!id>ShowerKinematics<!!id> object is attached to; the particle should be
  // the originator of the jet), and the momenta of the partners and 
  // the corresponding evolution scales for each type of interaction
  // considered (QCD, QED, EWK, ...) as defined in <!class>ShowerIndex<!!class>. 
  // This method should be used only if <!id>isTheJetStartingPoint()<!!id> is true.

private:

  ShowerKinematics & operator=(const ShowerKinematics &);
  //  Private and non-existent assignment operator.

  bool _isTheJetStartingPoint;

};

}

#include "ShowerKinematics.icc"

#endif /* HERWIG_ShowerKinematics_H */
