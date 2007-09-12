// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWHardProcess class.
//

#include "CKKWHardProcess.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWHardProcess.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CKKWHardProcess::~CKKWHardProcess() {}

void CKKWHardProcess::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void CKKWHardProcess::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

AbstractClassDescription<CKKWHardProcess> CKKWHardProcess::initCKKWHardProcess;
// Definition of the static class description member.

void CKKWHardProcess::Init() {

  static ClassDocumentation<CKKWHardProcess> documentation
    ("CKKWHardProcess is the base class for a boolean predicate"
     " specifying wether a hard process has been reached in reconstructing"
     " a parton shower history.");

}

