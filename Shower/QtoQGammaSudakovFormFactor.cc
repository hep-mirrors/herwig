// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSudakovFormFactor class.
//

#include "QtoQGammaSudakovFormFactor.h"

using namespace Herwig;


QtoQGammaSudakovFormFactor::~QtoQGammaSudakovFormFactor() {}


Energy QtoQGammaSudakovFormFactor::generateNextBranching( tPartCollHdlPtr ch, 
							  const Energy startingScale,
							  const bool reverseAngularOrder) {

  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to thie method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0; 

  //***LOOKHERE*** GENERATE  _q , AND EVENTUALLY ALSO  _z  AND  _phi

  return _q;

}

