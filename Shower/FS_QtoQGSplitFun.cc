// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtoQGSplitFun class.
//

#include "FS_QtoQGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


FS_QtoQGSplitFun::~FS_QtoQGSplitFun() {}


void FS_QtoQGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_QtoQGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_QtoQGSplitFun> FS_QtoQGSplitFun::initFS_QtoQGSplitFun;
// Definition of the static class description member.


void FS_QtoQGSplitFun::Init() {

  static ClassDocumentation<FS_QtoQGSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Final State NLO splitting function Q->QG.");

}

