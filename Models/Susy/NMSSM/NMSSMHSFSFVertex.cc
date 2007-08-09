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

void NMSSMHSFSFVertex::persistentOutput(PersistentOStream & os) const {
}

void NMSSMHSFSFVertex::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<NMSSMHSFSFVertex> NMSSMHSFSFVertex::initNMSSMHSFSFVertex;
// Definition of the static class description member.

void NMSSMHSFSFVertex::Init() {

  static ClassDocumentation<NMSSMHSFSFVertex> documentation
    ("There is no documentation for the NMSSMHSFSFVertex class");

}

void NMSSMHSFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				   tcPDPtr part2,tcPDPtr part3) {
}

void NMSSMHSFSFVertex::doinit() throw(InitException) {
  SSSVertex::doinit();
}
