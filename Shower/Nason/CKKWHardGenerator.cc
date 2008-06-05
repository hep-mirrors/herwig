// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWHardGenerator class.
//

#include "CKKWHardGenerator.h"
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
#include "Herwig++/Shower/Nason/NasonCKKWHandler.h"

using namespace Herwig;
using namespace ThePEG;

void CKKWHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _CKKWh;
}

void CKKWHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _CKKWh;
}

ClassDescription<CKKWHardGenerator> CKKWHardGenerator::initCKKWHardGenerator;
// Definition of the static class description member.

void CKKWHardGenerator::Init() {

  static ClassDocumentation<CKKWHardGenerator> documentation
    ("The CKKWHardGenerator class generates the hardest emission for"
     "vector boson decays to fermion-antifermion events in the Nason approach");

  static Reference<CKKWHardGenerator, NasonCKKWHandler> interfaceCKKWHandler
    ("CKKWHandler",
     "The cascade handler for NasonCKKW",
     &CKKWHardGenerator::_CKKWh, false, false,
     true, false);

}

NasonTreePtr CKKWHardGenerator::generateHardest(ShowerTreePtr tree) {
  //Get the NasonTree from the CKKW handler.
  NasonTreePtr nasontree = _CKKWh->getNasonTree();
	
  // Set the maximum pt for all other emissions
  // how should this be done for ckkw???
  // Energy ptveto(sqrt(_s)*_rt_mlambda/(4.*b_xs));
  // QProgenitor   ->maximumpT(ptveto);
  // QbarProgenitor->maximumpT(ptveto);

  // Connect the particles with the branchings in the NasonTree
  // nasontree->connect(QProgenitor->progenitor(),allBranchings[0]);
  // nasontree->connect(QbarProgenitor->progenitor(),allBranchings[1]);

  // Create the two colour lines and connect the particles:
  // ColinePtr blueLine  = new_ptr(ColourLine());
  // ColinePtr greenLine = new_ptr(ColourLine());
  // blueLine->addColoured(emitter,_iemitter);
  // blueLine->addColoured(gluon,_ispectator);
  // greenLine->addColoured(gluon,_iemitter);
  // greenLine->addColoured(parent,_iemitter);
  // greenLine->addColoured(spectator,_ispectator);

  // Calculate the shower variables:
  evolver()->showerModel()->kinematicsReconstructor()->
    reconstructDecayShower( nasontree, evolver() );

  // KMH - why don't we do the next step in reconstructDecayShower? 
  // Reset the momenta to ensure the correct momenta after shower recon
  // if emitter for Kleiss trick and shower are different
  for(map<ShowerParticlePtr,tNasonBranchingPtr>::const_iterator 
  	mit=nasontree->particles().begin();mit!=nasontree->particles().end();++mit)
    mit->first->set5Momentum(mit->second->_shower);

  // Return the NasonTree
  return nasontree;
}

//this is exactly as it was  for VectorBoson...
bool CKKWHardGenerator::canHandle(ShowerTreePtr tree) {
  if(tree->incomingLines().size()!=1) return false;    
  if((tree->incomingLines().begin()->first->id()==22)&&
  (tree->incomingLines().begin()->first->progenitor()->id()==23)) return false;

  //commented out the lines on outgoing particles since these are not currently
  //set correctly in the showerTree

  // map<ShowerProgenitorPtr,tShowerParticlePtr> outgoing=tree->outgoingLines();
  //  if(outgoing.size()!=2) return false;
  // if(abs(outgoing.begin()->first->progenitor()->id())>6)  return false;
  // if(outgoing.begin()->first->progenitor()->id()!=
  //   -1*outgoing.rbegin()->first->progenitor()->id())     return false;
  return true;
}

void CKKWHardGenerator::doinit() throw(InitException) {
  HardestEmissionGenerator::doinit();
}

void CKKWHardGenerator::dofinish() {
  HardestEmissionGenerator::dofinish();

}

void CKKWHardGenerator::doinitrun() {
  HardestEmissionGenerator::doinitrun();
}
