// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"

using namespace Herwig;


IS_QtildaShowerKinematics1to2::~IS_QtildaShowerKinematics1to2() {}


void IS_QtildaShowerKinematics1to2::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_QtildaShowerKinematics1to2::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_QtildaShowerKinematics1to2> 
IS_QtildaShowerKinematics1to2::initIS_QtildaShowerKinematics1to2;
// Definition of the static class description member.


void IS_QtildaShowerKinematics1to2::Init() {

  static ClassDocumentation<IS_QtildaShowerKinematics1to2> documentation
    ("This (concrete) class provides the specific Intial State shower kinematics information.");

}


void IS_QtildaShowerKinematics1to2::
updateChildren( const Energy qtilda, const double z, const double phi,
		const Energy Mass1, const Energy Mass2,
		double & sudAlpha1, double & px1, double &py1,
		double & sudAlpha2, double & px2, double &py2) {

  //***LOOKHERE*** WRITE THE CODE

}


void IS_QtildaShowerKinematics1to2::
updateParent( tCollecShoKinPtr & shoKinChildren ) {

  //***LOOKHERE*** WRITE THE CODE

}


Energy IS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}


