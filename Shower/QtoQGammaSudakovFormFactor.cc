// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSudakovFormFactor class.
//

#include "QtoQGammaSudakovFormFactor.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"

using namespace Herwig;


QtoQGammaSudakovFormFactor::~QtoQGammaSudakovFormFactor() {}


void QtoQGammaSudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void QtoQGammaSudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<QtoQGammaSudakovFormFactor> 
QtoQGammaSudakovFormFactor::initQtoQGammaSudakovFormFactor;
// Definition of the static class description member.


void QtoQGammaSudakovFormFactor::Init() {

  static ClassDocumentation<QtoQGammaSudakovFormFactor> documentation
    ("This (concrete) class describes the properties of ",
     "Sudakov form factor for Q->QGamma splitting.");

}


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

