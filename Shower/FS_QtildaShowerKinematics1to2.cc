// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtildaShowerKinematics1to2 class.
//

#include "FS_QtildaShowerKinematics1to2.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"

using namespace Herwig;


FS_QtildaShowerKinematics1to2::~FS_QtildaShowerKinematics1to2() {}


void FS_QtildaShowerKinematics1to2::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_QtildaShowerKinematics1to2::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_QtildaShowerKinematics1to2> FS_QtildaShowerKinematics1to2::initFS_QtildaShowerKinematics1to2;
// Definition of the static class description member.


void FS_QtildaShowerKinematics1to2::Init() {

  static ClassDocumentation<FS_QtildaShowerKinematics1to2> documentation
    ("This (concrete) class provides the specific Final State shower kinematics information.");

}


void FS_QtildaShowerKinematics1to2::
updateChildren( const Energy qtilda, const double z, const double phi,
		const Energy onShellMassChild1, const Energy onShellMassChild2,
		Ptr< QtildaShowerKinematics1to2 >::pointer & showerKinChild1,
		Ptr< QtildaShowerKinematics1to2 >::pointer & showerKinChild2 ) {

  //***LOOKHERE*** WRITE THE CODE
  showerKinChild1->alpha( alpha() * z );
  showerKinChild2->alpha( alpha() * (1.0 - z) );

  // Energy2 pPerp2 = sqr( 1.0 - z ) * ( sqr( z*qtilda ) - pVector().mass2() );

  // ... TO BE COMPLETED ...

}


void FS_QtildaShowerKinematics1to2::
updateParent( Ptr< QtildaShowerKinematics1to2 >::transient_pointer & showerKinChild1,
	      Ptr< QtildaShowerKinematics1to2 >::transient_pointer & showerKinChild2 ) {

  //***LOOKHERE*** WRITE THE CODE

}


Energy FS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}


