// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtildaShowerKinematics1to2 class.
//

#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"

using namespace Herwig;


FS_QtildaShowerKinematics1to2::~FS_QtildaShowerKinematics1to2() {}


void FS_QtildaShowerKinematics1to2::
updateChildren( const double parentSudAlpha, 
		const Energy parentSudPx, const Energy parentSudPy, 
		vector<double> & sudAlphaVect, 
		vector<Energy> & sudPxVect, vector<Energy> & sudPyVect ) {
 
  //***LOOKHERE*** WRITE THE CODE

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


