// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonEvolver class.
//

#include "NasonEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

void NasonEvolver::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void NasonEvolver::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<NasonEvolver> NasonEvolver::initNasonEvolver;
// Definition of the static class description member.

void NasonEvolver::Init() {

  static ClassDocumentation<NasonEvolver> documentation
    ("There is no documentation for the NasonEvolver class");

}

