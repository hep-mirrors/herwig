// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtoQGammaSplitFun class.
//

#include "FS_QtoQGammaSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


FS_QtoQGammaSplitFun::~FS_QtoQGammaSplitFun() {}


void FS_QtoQGammaSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void FS_QtoQGammaSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<FS_QtoQGammaSplitFun> FS_QtoQGammaSplitFun::initFS_QtoQGammaSplitFun;
// Definition of the static class description member.


void FS_QtoQGammaSplitFun::Init() {

  static ClassDocumentation<FS_QtoQGammaSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Final State NLO splitting function Q->QGamma.");

}

