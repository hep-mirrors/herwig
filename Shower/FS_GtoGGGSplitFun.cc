// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_GtoGGGSplitFun class.
//

#include "FS_GtoGGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


FS_GtoGGGSplitFun::~FS_GtoGGGSplitFun() {}


void FS_GtoGGGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_GtoGGGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_GtoGGGSplitFun> FS_GtoGGGSplitFun::initFS_GtoGGGSplitFun;
// Definition of the static class description member.


void FS_GtoGGGSplitFun::Init() {

  static ClassDocumentation<FS_GtoGGGSplitFun> documentation
    ("This is the (concrete) class which describes the Final State G->GGG splitting.");

}

