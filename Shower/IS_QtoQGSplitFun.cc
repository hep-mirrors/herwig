// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtoQGSplitFun class.
//

#include "IS_QtoQGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"


using namespace Herwig;


IS_QtoQGSplitFun::~IS_QtoQGSplitFun() {}


void IS_QtoQGSplitFun::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void IS_QtoQGSplitFun::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<IS_QtoQGSplitFun> IS_QtoQGSplitFun::initIS_QtoQGSplitFun;
// Definition of the static class description member.


void IS_QtoQGSplitFun::Init() {

  static ClassDocumentation<IS_QtoQGSplitFun> documentation
    ("This (concrete) class can (but it is not forced to) provide the exact ",
     "Initial State NLO splitting function Q->QG.");

}


