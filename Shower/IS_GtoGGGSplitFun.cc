// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_GtoGGGSplitFun class.
//

#include "IS_GtoGGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


IS_GtoGGGSplitFun::~IS_GtoGGGSplitFun() {}


void IS_GtoGGGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_GtoGGGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_GtoGGGSplitFun> IS_GtoGGGSplitFun::initIS_GtoGGGSplitFun;
// Definition of the static class description member.


void IS_GtoGGGSplitFun::Init() {

  static ClassDocumentation<IS_GtoGGGSplitFun> documentation
    ("This is the (concrete) class which describes the Initial State G->GGG splitting.");

}

