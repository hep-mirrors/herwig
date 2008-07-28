// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExternalHardGenerator class.
//

#include "ExternalHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

void ExternalHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _CKKWh;
}

void ExternalHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _CKKWh;
}

ClassDescription<ExternalHardGenerator> ExternalHardGenerator::initExternalHardGenerator;
// Definition of the static class description member.

void ExternalHardGenerator::Init() {

  static ClassDocumentation<ExternalHardGenerator> documentation
    ("The ExternalHardGenerator class generates the hardest emission for"
     "vector boson decays to fermion-antifermion events in the POWHEG approach");

  static Reference<ExternalHardGenerator, PowhegHandler> interfaceCKKWHandler
    ("CKKWHandler",
     "The cascade handler for Powheg",
     &ExternalHardGenerator::_CKKWh, false, false,
     true, false);

}

HardTreePtr ExternalHardGenerator::generateHardest(ShowerTreePtr tree) {
  
  ShowerProgenitorPtr 
    QProgenitor    = tree->outgoingLines().begin()->first,
    QbarProgenitor = tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor   ,QbarProgenitor);

  // Get the HardTree from the CKKW handler.
  HardTreePtr hardtree = _CKKWh->getHardTree();
  Energy kt_merge = _CKKWh->getMergeScale();

  //add pt vetos
  //n.b these are currently in terms of durham kt - should change veto
  QProgenitor->maximumpT(kt_merge);
  QbarProgenitor->maximumpT(kt_merge);

  if(!hardtree) return HardTreePtr();
  // check the trees match
  if(!hardtree->connect(tree)) return HardTreePtr();
  // Return the HardTree
  return hardtree;
}

bool ExternalHardGenerator::canHandle(ShowerTreePtr) {
  return true;
}
