// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSudakovFormFactor class.
//

#include "GtoGGSudakovFormFactor.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Repository/EventGenerator.h" 
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;


GtoGGSudakovFormFactor::~GtoGGSudakovFormFactor() {}


void GtoGGSudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void GtoGGSudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<GtoGGSudakovFormFactor> GtoGGSudakovFormFactor::initGtoGGSudakovFormFactor;
// Definition of the static class description member.


void GtoGGSudakovFormFactor::Init() {

  static ClassDocumentation<GtoGGSudakovFormFactor> documentation
    ("This (concrete) class describes the properties of ",
     "Sudakov form factor for G->GG splitting.");

}


Energy GtoGGSudakovFormFactor::generateNextBranching( tPartCollHdlPtr ch, 
						      const Energy startingScale,
						      const bool reverseAngularOrder) const {

  Energy newScale = Energy();

  //***LOOKHERE*** TO BE FILLED

  if (reverseAngularOrder) {
    newScale = startingScale / UseRandom::rnd();
  } else {
    newScale = startingScale * UseRandom::rnd(); 
  }
  
  // first toy

  return newScale;

}

