// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_GtoQQbarSplitFun class.
//

#include "IS_GtoQQbarSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


IS_GtoQQbarSplitFun::~IS_GtoQQbarSplitFun() {}


void IS_GtoQQbarSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_GtoQQbarSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_GtoQQbarSplitFun> IS_GtoQQbarSplitFun::initIS_GtoQQbarSplitFun;
// Definition of the static class description member.


void IS_GtoQQbarSplitFun::Init() {

  static ClassDocumentation<IS_GtoQQbarSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Initial State NLO splitting function G->QQbar.");

}

