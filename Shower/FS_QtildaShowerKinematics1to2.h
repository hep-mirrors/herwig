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
					const Lorentz5Momentum & n,
					const Energy inputOnShellMass );
  // Creator with the two defining vectors <!id>p<!!id> and <!id>n<!!id> 
  // and the on-shell mass of the particle.

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
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> 
  // is true, and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods

private:

  static ClassDescription<FS_QtildaShowerKinematics1to2> initFS_QtildaShowerKinematics1to2;
  // Describe an abstract base class with persistent data.

  FS_QtildaShowerKinematics1to2 & operator=(const FS_QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of FS_QtildaShowerKinematics1to2.
template <>
struct BaseClassTrait<Herwig::FS_QtildaShowerKinematics1to2,1> {
  typedef Herwig::QtildaShowerKinematics1to2 NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::FS_QtildaShowerKinematics1to2>: public ClassTraitsBase<Herwig::FS_QtildaShowerKinematics1to2> {
  static string className() { return "/Herwig++/FS_QtildaShowerKinematics1to2"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "FS_QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_FS_QtildaShowerKinematics1to2_H */
