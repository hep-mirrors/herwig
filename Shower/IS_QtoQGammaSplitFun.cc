// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtoQGammaSplitFun class.
//

#include "IS_QtoQGammaSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


IS_QtoQGammaSplitFun::~IS_QtoQGammaSplitFun() {}


void IS_QtoQGammaSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_QtoQGammaSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_QtoQGammaSplitFun> IS_QtoQGammaSplitFun::initIS_QtoQGammaSplitFun;
// Definition of the static class description member.


void IS_QtoQGammaSplitFun::Init() {

  static ClassDocumentation<IS_QtoQGammaSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Initial State NLO splitting function Q->QGamma.");

}

