// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtildaShowerKinematics1to2 class.
//

#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Repository/EventGenerator.h"

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
		const Energy Mass1, const Energy Mass2,
		double & sudAlpha1, double & px1, double & py1,
		double & sudAlpha2, double & px2, double & py2) {
  
  //***LOOKHERE*** complete the code  
  sudAlpha1 *= z; 
  sudAlpha2 *= 1.-z; 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl; 
    if ( sqr( z*qtilda ) - sqr(Mass1) < 0 ) {
      generator()->log() << "Warning! check phase space of qtilda and z!" << endl; 
    }
  }

  Energy pPerp = (1.-z)*sqrt( sqr( z*qtilda ) - sqr(Mass1));
  px1 = pPerp*cos(phi) + z*px1; 
  py1 = pPerp*sin(phi) + z*py1; 
  px2 -= px1; 
  py2 -= py1;

  // ... TO BE COMPLETED AND checked!...

}


void FS_QtildaShowerKinematics1to2::
updateParent( tCollecShoKinPtr & shoKinChildren ) {

  //***LOOKHERE*** WRITE THE CODE

}


Energy FS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}


