// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSudakovFormFactor class.
//

#include "GtoGGSudakovFormFactor.h"
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;


GtoGGSudakovFormFactor::~GtoGGSudakovFormFactor() {}


Energy GtoGGSudakovFormFactor::generateNextBranching( tPartCollHdlPtr ch, 
						      const Energy startingScale,
						      const bool reverseAngularOrder) {

  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to thie method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0; 

  //***LOOKHERE*** GENERATE  _q , AND EVENTUALLY ALSO  _z  AND  _phi
  //               BELOW IS JUST A TEMPORARY FAKE
  if (reverseAngularOrder) {
    _q = startingScale / UseRandom::rnd();
  } else {
    _q = startingScale * UseRandom::rnd();
  }
  _z = 0.0;
  _phi = ( UseRandom::rndbool() ? 1.0 : -1.0 ) * 3.1415*UseRandom::rnd();
 
  return _q;

}

