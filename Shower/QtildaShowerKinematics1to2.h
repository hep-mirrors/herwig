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
//   *** LOOKHERE***  VERY PRELIMINARY!!! 
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

  virtual void updateChildren( const Energy qtilda, const double z, const double phi,
			       const Energy onShellMassChild1, const Energy onShellMassChild2,
			       Ptr< QtildaShowerKinematics1to2 >::pointer & showerKinChild1,
			       Ptr< QtildaShowerKinematics1to2 >::pointer & showerKinChild2 ) = 0;
  // Along with the showering evolution --- 
  // going forward for time-like (forward) evolution, and
  // going backward for space-like (backward) evolution ---
  // the shower kinematics of the branching products is created and filled.

  virtual void updateParent( Ptr< QtildaShowerKinematics1to2 >::transient_pointer & showerKinChild1,
			     Ptr< QtildaShowerKinematics1to2 >::transient_pointer & showerKinChild2 ) = 0;
  // After showering, moving in the reverse direction ---
  // going backward for time-like (forward) evolution, and
  // going forward for space-like (backward) evolution ---
  // the shower kinematics of the parent is updated (completed)
  // from the kinematics of its children.

  virtual Energy jetMass() = 0;
  // Pure virtual method, to be defined in a derived class.
  // The method returns the mass of jet. 
  // It should be used only if <!id>isTheJetStartingPoint()<!!id> is true, 
  // and only after the jet kinematics reconstruction.
  // (performed by the <!class>KinematicsReconstructor<!!class> class).

  inline Energy qtilde() const;
  inline void qtilde( const Energy inputQtilde );
  inline double z() const;
  inline void z( const double inputZ );
  inline double phi() const;
  inline void phi( const double inputPhi );
  // Access/set to the generated kinematics variables of the splitting <I>1-&GT;2</I>.

  inline double alpha() const;
  inline void alpha( const double inputAlpha );
  inline double beta() const;
  inline void beta( const double inputBeta );
  inline Energy qPerp1() const;
  inline void qPerp1( const Energy inputQPerp1 );
  inline Energy qPerp2() const;
  inline void qPerp2( const Energy inputQPerp2 );
  inline Energy pPerp1() const;
  inline void pPerp1( const Energy inputPPerp1 );
  inline Energy pPerp2() const;
  inline void pPerp2( const Energy inputPPerp2 );
  // Access/set to the calculated kinematics variables of the splitting <I>1-&GT;2</I>.

  inline Energy onShellMass() const;
  inline void onShellMass( const Energy inputOnShellMass );

  inline const Lorentz5Momentum & pVector() const;
  inline const Lorentz5Momentum & nVector() const;
  // Access to the <!id>p<!!id> and <!id>n<!!id> vectors used to described the kinematics.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static AbstractClassDescription<QtildaShowerKinematics1to2> initQtildaShowerKinematics1to2;
  // Describe an abstract base class with persistent data.

  QtildaShowerKinematics1to2 & operator=(const QtildaShowerKinematics1to2 &);
  //  Private and non-existent assignment operator.

  Energy _qtilde;
  double _z;
  double _phi;
  double _alpha;
  double _beta;
  Energy _qPerp1;
  Energy _qPerp2;
  Energy _pPerp1;
  Energy _pPerp2;
  Energy _onShellMass;
  Lorentz5Momentum _pVector;
  Lorentz5Momentum _nVector;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of QtildaShowerKinematics1to2.
template <>
struct BaseClassTrait<Herwig::QtildaShowerKinematics1to2,1> {
  typedef Herwig::ShowerKinematics NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::QtildaShowerKinematics1to2>: public ClassTraitsBase<Herwig::QtildaShowerKinematics1to2> {
  static string className() { return "/Herwig++/QtildaShowerKinematics1to2"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_QtildaShowerKinematics1to2_H */
