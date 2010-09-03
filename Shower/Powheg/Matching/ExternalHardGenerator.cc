// -*- C++ -*-//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExternalHardGenerator class.
//

#include "ExternalHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
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
#include "Herwig++/Shower/Base/CKKWVeto.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

void ExternalHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << CKKWh_;
}

void ExternalHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> CKKWh_;
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
     &ExternalHardGenerator::CKKWh_, false, false,
     true, false);

}

HardTreePtr ExternalHardGenerator::generateHardest(ShowerTreePtr tree) {
  //get com energy from the progenitors
  Energy2 s = 
    (tree->incomingLines(). begin()->first->progenitor()->momentum() +
     tree->incomingLines().rbegin()->first->progenitor()->momentum()).m2();
  // Get the HardTree from the CKKW handler.
  CKKWTreePtr hardtree = CKKWh_->getCKKWTree();
  Energy veto_pt = sqrt( s );
  // evolver vetoes should be set to use shower pt in highest multiplicity case
  if( CKKWh_->highestMult() && hardtree ) veto_pt = hardtree->lowestPt( 1, s );
  for( map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it
	 = tree->outgoingLines().begin(); 
       it != tree->outgoingLines().end(); ++it ) {
    if( ! it->second->coloured() ) continue;
    it->first->maximumpT( veto_pt );
  }     
  for( map< ShowerProgenitorPtr, ShowerParticlePtr >::iterator it
	 = tree->incomingLines().begin(); 
       it != tree->incomingLines().end(); ++it ) {
    if( ! it->second->coloured() ) continue;
    it->first->maximumpT( veto_pt );
  }   
  // if no hard tree
  if( ! hardtree ) return HardTreePtr();
  // check the trees match
  if( ! hardtree->connect(tree) ) return HardTreePtr();
  // if not lowest multiplcity return
  if( !CKKWh_->lowestMult() ) return hardtree; 
  // if no correction return
  if( !CKKWh_->matrixElement() || 
      !CKKWh_->matrixElement()->hasPOWHEGCorrection() ) return hardtree;
  // fill the dead zone by generating an emission with the hard generator
  // updating the _theHardTree and vetoeing the event if the emission pt 
  // (from lowest pt) is greater than the merge scale
  HardTreePtr newHardTree = CKKWh_->matrixElement()->generateHardest( tree );
  if(!newHardTree) return hardtree;
  ShowerParticleVector particles;
  // Sudakovs for ISR
  for(set<HardBranchingPtr>::iterator cit=newHardTree->branchings().begin();
      cit!=newHardTree->branchings().end();++cit) {
    if((**cit).parent()&&(**cit).status()==HardBranching::Incoming) {
      IdList br(3);
      br[0] = (**cit).parent()->branchingParticle()->id();
      br[1] = (**cit).          branchingParticle()->id();
      br[2] = (**cit).parent()->children()[0]==*cit ?
	(**cit).parent()->children()[1]->branchingParticle()->id() :
	(**cit).parent()->children()[0]->branchingParticle()->id();
      BranchingList branchings = 
	evolver()->splittingGenerator()->initialStateBranchings();
      if(br[1]<0&&br[0]==br[1]) {
	br[0] = abs(br[0]);
	br[1] = abs(br[1]);
      }
      else if(br[1]<0) {
	br[1] = -br[1];
	br[2] = -br[2];
      }
      long index = abs(br[1]);
      SudakovPtr sudakov;
      for(BranchingList::const_iterator cjt = branchings.lower_bound(index); 
	  cjt != branchings.upper_bound(index); ++cjt ) {
	IdList ids = cjt->second.second;
	if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
	  sudakov=cjt->second.first;
	  break;
	}
      }
      if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				     << "Evolver::generateHardest()" 
				     << Exception::runerror;
      (**cit).parent()->sudakov(sudakov);
    }
    else if(!(**cit).children().empty()) {
      IdList br(3);
      br[0] = (**cit)               .branchingParticle()->id();
      br[1] = (**cit).children()[0]->branchingParticle()->id();
      br[2] = (**cit).children()[1]->branchingParticle()->id();
      BranchingList branchings = 
	evolver()->splittingGenerator()->finalStateBranchings();
      if(br[0]<0) {
	br[0] = abs(br[0]);
	br[1] = abs(br[1]);
	br[2] = abs(br[2]);
      }
      long index = br[0];
      SudakovPtr sudakov;
      for(BranchingList::const_iterator cjt = branchings.lower_bound(index); 
	  cjt != branchings.upper_bound(index); ++cjt ) {
	IdList ids = cjt->second.second;
	if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
	  sudakov=cjt->second.first;
	  break;
	}
      }
      if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				     << "Evolver::generateHardest()" 
				     << Exception::runerror;
      (**cit).sudakov(sudakov);
    }
  }
  // calculate the evolution scale
  for(set<HardBranchingPtr>::iterator cit=newHardTree->branchings().begin();
      cit!=newHardTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  evolver()->showerModel()->
    partnerFinder()->setInitialEvolutionScales(particles,false,
					       newHardTree->interaction(),true);
  // inverse reconstruction
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructHardJets(newHardTree,ShowerHandler::currentHandler()->evolver(),
			newHardTree->interaction());

  // \todo needs sorting out !!!!!!!!!!

//   newHardTree->findNodes();
//   // reject the event if newnewHardTree->lowestPt
//   // is greater than pt_cut
//   if( newHardTree->lowestPt( 1, s ) > CKKWh_->theVeto()->getVeto() ) 
//     throw Veto();



  //maximum pt of shower progenitors will also have been 
  //updated to something below the veto scale
  //set the CKKWVeto to use the maximum pt scale as veto scale 
  //rather than the fixed veto scale
  //this means doing it as if it is a highest multiplicity emission 
  CKKWh_->theVeto()->setHighest( true );
  return newHardTree;
}

bool ExternalHardGenerator::canHandle(ShowerTreePtr) {
  return true;
}
