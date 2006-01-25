// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerKinematics class.
//

#include "ShowerKinematics.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ShowerIndex.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerKinematics.tcc"
#endif


using namespace Herwig;

ShowerKinematics::~ShowerKinematics() {}

Lorentz5Momentum ShowerKinematics::
referenceFrame( const Lorentz5Momentum & particleMomentum,
		const vector<Lorentz5Momentum> & partnersMomenta,
		const vector<Energy> & evolutionScales ) {

  Lorentz5Momentum cmFrame = Lorentz5Momentum();

  //***LOOKHERE*** WRITE THE CODE!

  return cmFrame;

}

