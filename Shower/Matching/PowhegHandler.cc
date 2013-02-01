// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegHandler class.
//

#include "PowhegHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/PDF/PDF.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"

using namespace Herwig;

IBPtr PowhegHandler::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegHandler::fullclone() const {
  return new_ptr(*this);
}

void PowhegHandler::persistentOutput(PersistentOStream & os) const {
  os << pTDefinition_ << ounit(maxpT_,GeV);
}

void PowhegHandler::persistentInput(PersistentIStream & is, int) {
  is >> pTDefinition_ >> iunit(maxpT_,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegHandler,MatchingHandler>
describeHerwigPowhegHandler("Herwig::PowhegHandler", "HwMatching.so");

void PowhegHandler::Init() {

  static ClassDocumentation<PowhegHandler> documentation
    ("The PowhegHandler class handles the generation of the shower, including truncated"
     "showers for external processes using the POWHEG scheme.");

  static Switch<PowhegHandler,unsigned int> interfacepTDefinition
    ("pTDefinition",
     "The choice of the definition of the pT for the maximum emission in the shower",
     &PowhegHandler::pTDefinition_, 0, false, false);
  static SwitchOption interfacepTDefinitionScale
    (interfacepTDefinition,
     "Scale",
     "Use the value of the lastScale(), which if reading Les Houches events is SCALUP",
     0);
  static SwitchOption interfacepTDefinitionScaleAndpT
    (interfacepTDefinition,
     "ScaleAndpT",
     "Use the same choice as scale for non-emission events but the pT of the"
     " emitted particle for emission events",
     1);
  static SwitchOption interfacepTDefinitionmaxpTAndpT
    (interfacepTDefinition,
     "maxpTAndpT",
     "Use the maxPt value fo non-emission events and the pT of the"
     " emitted particle for emission events",
     2);

  static Parameter<PowhegHandler,Energy> interfacemaxpT
    ("maxpT",
     "Maximum pT for emission from non-emission events.",
     &PowhegHandler::maxpT_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

double PowhegHandler::reweightCKKW(int minMult, int maxMult) {
  // as not reweighting skip if in initialisation
  if(generator()->state()==InterfacedBase::initializing)
    return 1.;
  initialiseMatching(minMult,maxMult);
  if( ! hardTree().tree() ) return 1.;
  // update sub process based on hardTree()
  updateSubProcess();
  return 1.;
}

HardTreePtr PowhegHandler::generateCKKW(ShowerTreePtr tree) const {
  //get com energy from the progenitors
  Energy2 s = 
    (tree->incomingLines(). begin()->first->progenitor()->momentum() +
     tree->incomingLines().rbegin()->first->progenitor()->momentum()).m2();
  // Get the HardTree from the CKKW handler.
  CKKWTreePtr hardtree = getCKKWTree();
  // if no hard tree return
  if( ! hardtree ) return HardTreePtr();
  // check the trees match
  if( ! hardtree->connect(tree) ) return HardTreePtr();
  // choice of the veto pT
  Energy veto_pt = sqrt( s );
  // use the scale
  if(pTDefinition_==0) {
    veto_pt = sqrt(lastXCombPtr()->lastScale());
  }
  // otherwise pT of emission if emitted
  else if(!lowestMult()) {
    veto_pt = hardtree->lowestPt( 1, s );
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
  }
  else if(pTDefinition_==1) {
    veto_pt = sqrt(lastXCombPtr()->lastScale());
  }
  else if(pTDefinition_==2) {
    veto_pt = maxpT_;
  }
  else
    assert(false);
  // that's it
  return hardtree;
}

PotentialTree PowhegHandler::chooseHardTree(double,double) {
  PotentialTree chosen_hardTree;
  Energy min_pt = Constants::MaxEnergy;
  if( !hardTrees().empty() ){
    for(unsigned int ix = 0; ix < hardTrees().size(); ix++ ){
      if( hardTrees()[ix].first.tree()->totalPt() < min_pt ) {
	min_pt = hardTrees()[ix].first.tree()->totalPt();
	chosen_hardTree = hardTrees()[ix].first;
      }
    }
  }
  else {
    for(unsigned int ix = 0; ix < nonOrderedTrees().size(); ix++ ){
      if( nonOrderedTrees()[ix].first.tree()->totalPt() < min_pt ){
	min_pt = nonOrderedTrees()[ix].first.tree()->totalPt();
	chosen_hardTree=nonOrderedTrees()[ix].first;
      }
    }
  }
  return chosen_hardTree;
}
