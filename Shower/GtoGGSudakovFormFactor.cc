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
    _z = UseRandom::rnd(); 
  } else {
    get_qz(true, -1.1, .15, max(750.*MeV, splitFun()->massEmitter()), startingScale, _q, _z); 
  }
  _phi = 2.*pi*UseRandom::rnd();
 
  return _q;

}

