// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerKinematics class.
//

#include "ShowerKinematics.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/CLHEPWrap/Lorentz5Vector.h"
#include "ShowerIndex.h"

using namespace Herwig;


ShowerKinematics::~ShowerKinematics() {}


AbstractClassDescription<ShowerKinematics> ShowerKinematics::initShowerKinematics;
// Definition of the static class description member.


void ShowerKinematics::Init() {

  static ClassDocumentation<ShowerKinematics> documentation
    ("This abstract class is the base for all other shower kinematics classes.");

}


Lorentz5Momentum ShowerKinematics::
referenceFrame( const Lorentz5Momentum & particleMomentum,
		const vector<Lorentz5Momentum> & partnersMomenta,
		const vector<Energy> & evolutionScales ) {

  Lorentz5Momentum cmFrame = Lorentz5Momentum();

  //***LOOKHERE*** WRITE THE CODE!

  return cmFrame;

}



