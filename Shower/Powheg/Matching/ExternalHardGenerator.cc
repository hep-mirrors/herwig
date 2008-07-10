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
  //Get the HardTree from the CKKW handler.
  HardTreePtr hardtree = _CKKWh->getHardTree();
  vector<ShowerProgenitorPtr> progenitors = tree->extractProgenitors();
  // connect the trees up
  for(unsigned int ix=0;ix<progenitors.size();++ix) {
    bool match(false);
    for(set<HardBranchingPtr>::iterator it=hardtree->branchings().begin();
	it!=hardtree->branchings().end();++it) {
      if((**it).branchingParticle()->id()!=progenitors[ix]->progenitor()->id()) continue;
      cerr << *progenitors[ix]->progenitor() << "\n";
      hardtree->connect(progenitors[ix]->progenitor(),*it);
      match = true;
    }
    if(!match) return HardTreePtr();
  }











  // Set the maximum pt for all other emissions
  // how should this be done for ckkw???
  // Energy ptveto(sqrt(_s)*_rt_mlambda/(4.*b_xs));
  // QProgenitor   ->maximumpT(ptveto);
  // QbarProgenitor->maximumpT(ptveto);

  // Connect the particles with the branchings in the HardTree
  // hardtree->connect(QProgenitor->progenitor(),allBranchings[0]);
  // hardtree->connect(QbarProgenitor->progenitor(),allBranchings[1]);

  // Create the two colour lines and connect the particles:
  // ColinePtr blueLine  = new_ptr(ColourLine());
  // ColinePtr greenLine = new_ptr(ColourLine());
  // blueLine->addColoured(emitter,_iemitter);
  // blueLine->addColoured(gluon,_ispectator);
  // greenLine->addColoured(gluon,_iemitter);
  // greenLine->addColoured(parent,_iemitter);
  // greenLine->addColoured(spectator,_ispectator);

  // Return the HardTree
  return hardtree;
}

bool ExternalHardGenerator::canHandle(ShowerTreePtr tree) {
  return true;
}
