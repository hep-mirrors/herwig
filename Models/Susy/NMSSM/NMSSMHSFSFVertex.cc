// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHSFSFVertex class.
//

#include "NMSSMHSFSFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMHSFSFVertex::NMSSMHSFSFVertex() {}

void NMSSMHSFSFVertex::persistentOutput(PersistentOStream & ) const {
}

void NMSSMHSFSFVertex::persistentInput(PersistentIStream & , int) {
}

ClassDescription<NMSSMHSFSFVertex> NMSSMHSFSFVertex::initNMSSMHSFSFVertex;
// Definition of the static class description member.

void NMSSMHSFSFVertex::Init() {

  static ClassDocumentation<NMSSMHSFSFVertex> documentation
    ("There is no documentation for the NMSSMHSFSFVertex class");

}

void NMSSMHSFSFVertex::setCoupling(Energy2 ,tcPDPtr ,
				   tcPDPtr,tcPDPtr) {
}

void NMSSMHSFSFVertex::doinit() throw(InitException) {
  SSSVertex::doinit();
}
