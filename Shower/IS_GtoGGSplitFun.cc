// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_GtoGGSplitFun class.
//

#include "IS_GtoGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


IS_GtoGGSplitFun::~IS_GtoGGSplitFun() {}


void IS_GtoGGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_GtoGGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_GtoGGSplitFun> IS_GtoGGSplitFun::initIS_GtoGGSplitFun;
// Definition of the static class description member.


void IS_GtoGGSplitFun::Init() {

  static ClassDocumentation<IS_GtoGGSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Initial State NLO splitting function G->GG.");

}

