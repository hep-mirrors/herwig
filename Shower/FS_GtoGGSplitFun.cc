// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_GtoGGSplitFun class.
//

#include "FS_GtoGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


FS_GtoGGSplitFun::~FS_GtoGGSplitFun() {}


void FS_GtoGGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_GtoGGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_GtoGGSplitFun> FS_GtoGGSplitFun::initFS_GtoGGSplitFun;
// Definition of the static class description member.


void FS_GtoGGSplitFun::Init() {

  static ClassDocumentation<FS_GtoGGSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Final State NLO splitting function G->GG.");

}

