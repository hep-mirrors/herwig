// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_GtoQQbarSplitFun class.
//

#include "FS_GtoQQbarSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


FS_GtoQQbarSplitFun::~FS_GtoQQbarSplitFun() {}


void FS_GtoQQbarSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_GtoQQbarSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_GtoQQbarSplitFun> FS_GtoQQbarSplitFun::initFS_GtoQQbarSplitFun;
// Definition of the static class description member.


void FS_GtoQQbarSplitFun::Init() {

  static ClassDocumentation<FS_GtoQQbarSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Final State NLO splitting function G->QQbar.");

}

