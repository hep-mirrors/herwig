// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSudakovFormFactor class.
//

#include "QtoQGSudakovFormFactor.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Repository/EventGenerator.h" 
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;


QtoQGSudakovFormFactor::~QtoQGSudakovFormFactor() {}


void QtoQGSudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void QtoQGSudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<QtoQGSudakovFormFactor> QtoQGSudakovFormFactor::initQtoQGSudakovFormFactor;
// Definition of the static class description member.


void QtoQGSudakovFormFactor::Init() {

  static ClassDocumentation<QtoQGSudakovFormFactor> documentation
    ("This (concrete) class describes the properties of ",
     "Sudakov form factor for Q->QG splitting.");

}


Energy QtoQGSudakovFormFactor::generateNextBranching( tPartCollHdlPtr ch, 
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

