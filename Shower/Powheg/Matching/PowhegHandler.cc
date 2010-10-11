#include "PowhegHandler.h"
#include "ThePEG/Utilities/CFileLineReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include <queue>
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/PDF/PDF.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

IBPtr PowhegHandler::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegHandler::fullclone() const {
  return new_ptr(*this);
}

void PowhegHandler::persistentOutput(PersistentOStream & os) const {
  os  << _alphaS << _jetMeasureMode << _allowedInitial
      << _allowedFinal << _matrixElement << _highestMultOption 
      << _reweightOpt << _clusterOption << HWmatrixElement_
      << _rejectNonAO << _rejectNoHistories 
      << partonExtractor_ << cuts_ << _showerVeto << _cutOption;
}

void PowhegHandler::persistentInput(PersistentIStream & is, int) {
  is  >> _alphaS >> _jetMeasureMode >> _allowedInitial
      >> _allowedFinal >> _matrixElement >> _highestMultOption 
      >> _reweightOpt >> _clusterOption >> HWmatrixElement_ 
      >> _rejectNonAO >> _rejectNoHistories 
      >> partonExtractor_ >> cuts_ >> _showerVeto >> _cutOption;
}

ClassDescription<PowhegHandler> PowhegHandler::initPowhegHandler;
// Definition of the static class description member.

void PowhegHandler::Init() {

  static ClassDocumentation<PowhegHandler> documentation
    ("The PowhegHandler class manages the implementation of the CKKW approach using"
     "the truncated shower.");

  static Reference<PowhegHandler,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &PowhegHandler::_alphaS, false, false, true, false, false);

  static Reference<PowhegHandler,PartonExtractor> interfacePartonExtractor
    ("PartonExtractor",
     "The PartonExtractor object used to construct remnants. If no object is "
     "provided the LesHouchesEventHandler object must provide one instead.",
     &PowhegHandler::partonExtractor_, true, false, true, true, false);

  static Switch<PowhegHandler, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &PowhegHandler::_jetMeasureMode, 1, false, false);        
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  
  static SwitchOption ShowerPt
    (ifaceJetMeasureMode,"ShowerPt","ShowerPt", 1);
  
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 2);

  static Switch<PowhegHandler,bool> interfaceHighestMultiplicity
     ("HighestMultOption",
     "Whether we are using the highest mutliplicity treatment",
     &PowhegHandler::_highestMultOption, true, false, false);
  static SwitchOption interfaceMultHighest
    (interfaceHighestMultiplicity,
     "Yes",
     "Use the special option",
     true);
  static SwitchOption interfaceMultNotHighest
    (interfaceHighestMultiplicity,
     "No",
     "Don'y use the special option",
     false);

  static Switch<PowhegHandler, unsigned int> interfaceReweight
    ("ReweightOption",
     "Whether to switch off the sudakov reweighting",
     &PowhegHandler::_reweightOpt, 0, false, false);
  static SwitchOption interfaceReweightOff
    (interfaceReweight,
     "Off",
     "No Sudakov reweighting",
     1);
  static SwitchOption interfaceReweightNoSud
    (interfaceReweight,
     "NoSud",
     "alphaS and pdf reweighting only",
     2);

  static Reference<PowhegHandler,MEBase> interfaceMatrixElement
    ("MatrixElement",
     "The matrix element class for the core 2->2 process",
     &PowhegHandler::_matrixElement, false, false, true, false, false);

  static Switch<PowhegHandler, unsigned int> ifaceClusterOption
    ("ClusterOption",
     "Choice of the clustering scheme",
     &PowhegHandler::_clusterOption, 0, false, false);
  static SwitchOption allHistories
    (ifaceClusterOption,"allHistories","make all histories and require angular ordering", 0);
  static SwitchOption jetClustered
    (ifaceClusterOption,"ptChoice", "choose ordered history with lowest total pt", 2);
  static SwitchOption highProbChoice
    (ifaceClusterOption,"highProbChoice", "choose ordered history with highest probability", 3);

 static Switch<PowhegHandler,bool> interfaceReject
    ("RejectNonOrdered",
     "Whether to reject non angular ordered cluster histories",
     &PowhegHandler::_rejectNonAO, true, false, false);
 static SwitchOption interfaceRejectYes
    (interfaceReject,
     "Reject",
     "Reject all non angular ordered events",
     true);
  static SwitchOption interfaceRejectNo
    (interfaceReject,
     "Select",
     "Select a history for non angular ordered events",
     false);

  static Switch<PowhegHandler,bool> interfaceRejectNoHist
    ("RejectNoHistories",
     "Whether to reject events with no shower interpretation",
     &PowhegHandler::_rejectNoHistories, true, false, false);
  static SwitchOption interfaceRejectNoHistYes
    (interfaceRejectNoHist,
     "Reject",
     "Reject all events with no shower interpretation",
     true);
  static SwitchOption interfaceRejectNoHistNo
    (interfaceRejectNoHist,
     "Shower",
     "Shower events with no shower interpretation directly",
     false);

  static Reference<PowhegHandler,Cuts> interfaceCuts
    ("Cuts",
     "The Cuts object to be used for this reader. Note that these "
     "must not be looser cuts than those used in the actual generation. "
     "If no object is provided the LesHouchesEventHandler object must "
     "provide one instead.",
     &PowhegHandler::cuts_, true, false, true, true, false);
 
   static Reference<PowhegHandler,CKKWVeto > interfaceShowerVeto
    ("ShowerVeto",
     "The veto to be applied to the shower",
     &PowhegHandler::_showerVeto, false, false, true, false ); 

   static Switch<PowhegHandler, int> ifaceDynamicCutOption
     ("DynamicCutOption",
      "Options controlling whether dynamic cuts are applied to"
      " ensure that the matrix element configuration is resolved at the merging scale",
     &PowhegHandler::_cutOption, 0, false, false);
  static SwitchOption noDynamicCut
    (ifaceDynamicCutOption,"noCut","no dynamic cuts applied", 0);
  static SwitchOption allDynamicCut
    (ifaceDynamicCutOption,"cutAll", "apply resolution cut to full configuration", 1);
  static SwitchOption externalDynamicCut
    (ifaceDynamicCutOption,"cutExternal", "apply resolution cut to external partons only", 2);

}

double PowhegHandler::reweightCKKW(int minMult, int maxMult) {
  // check multiplicity of FS partons
  int nOut = lastXCombPtr()->subProcess()->outgoing().size();
  _highestMult = nOut == maxMult && _highestMultOption;
  _lowestMult  = nOut == minMult;
  _s = lastXCombPtr()->lastSHat();
  _pdfScale = sqrt(lastXCombPtr()->lastScale());
  _global_alphaS_wgt = 10.;
  _alphaSMG = lastXCombPtr()->lastAlphaS();
  // create a hard tree by clustering the event
  hardTree_ = PotentialTree();
  hardTree_ = doClustering();
  // set whether this is the highest multiplicity channel
  _showerVeto->setHighest( _highestMult );
  // set the veto to veto events or emissions depending
  // on whether Sudakovs are being generated dynamically
  _showerVeto->setDynamicSuds( ! _noShowerHists );
  // if no shower history is found either throw away or shower directly
  if( ! hardTree_.tree ) return _noShowerHists ? 1. : 0.;
  // apply dynamic matrix element vetoes
  if( _cutOption != 0 && 
      hardTree_.tree->lowestPtMomentum( _jetMeasureMode, _cutOption )  < _showerVeto->getVeto() ) 
    return 0.;
  // update sub process based on hardTree_
  updateSubProcess();
  // compute the Sudakov and alphaS weight
  if( _reweightOpt == 1) return 1.;
  hardTree_.tree->findNodes();
  double SudWgt = sudakovWeight( hardTree_.tree );    
  assert( !isnan( SudWgt ) && !isinf( SudWgt ) );
  return SudWgt;
}

void PowhegHandler::dofinish() {
  ShowerHandler::dofinish();
  if( _clusterOption == 2 ){
    generator()->log() <<"\n-----Powheg Handler----\n"
		       <<"proportion of events with no ordered trees created = "
		       << double( _unorderedEvents / double( _trees_created ) ) * 100.
		       <<" %\n---------------------\n";
  }
}

void PowhegHandler::doinit() {
  ShowerHandler::doinit();
  // extract the allowed branchings
  // final-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->finalStateBranchings().begin();
      it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    pair<long,long> prod(make_pair(it->second.second[1],it->second.second[2]));
    _allowedFinal.insert(make_pair(prod,it->second));
    swap(prod.first,prod.second);
    _allowedFinal.insert(make_pair(prod,it->second));
  }
  // initial-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->initialStateBranchings().begin();
      it != evolver()->splittingGenerator()->initialStateBranchings().end(); ++it) {
    _allowedInitial.insert(make_pair(it->second.second[0],it->second));
  }
  HWmatrixElement_ = dynamic_ptr_cast<HwMEBasePtr>(_matrixElement);
}

//given two particles returns value of durham jet algorithm
bool PowhegHandler::splittingAllowed( ShowerParticlePtr part_i,
				      ShowerParticlePtr part_j ) {
  //only consider QCD clusterings
  if( !( abs ( part_i->id() ) < 7 || part_i->id() == 21 ) ||
      !( abs ( part_j->id() ) < 7 || part_j->id() == 21 ) ) return false;
  //qqbar clustering checks
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() ) ) return false;
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return false;
  }
  return true;
}

// finds the id of the emitting particle and sudakov for the desired clustering
// also swaps order of children pointers as required (not pointers passed by reference)
SudakovPtr PowhegHandler::getSud( long & emmitter_id,
				  ShowerParticlePtr & part_i, 
				  ShowerParticlePtr & part_j ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() )  ) return SudakovPtr();
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return SudakovPtr();
    //if the q and qbar are the wrong way round then switch order
    if ( part_j->id() > part_i->id() ) swap( part_i, part_j );
    emmitter_id = 21;
  }
  // q/qbar 2 q/qbar g
  else if ( abs ( part_i->id() ) < 7 || abs ( part_j->id() ) < 7 ) {
    if( abs ( part_i->id() ) < 7 ){
      emmitter_id = part_i->id();
    }
    else {
      emmitter_id = part_j->id();
    }
  }
  // g 2 g g
  else {
    emmitter_id = 21;
  }
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();

  //cycle through list of all branchings with the correct abs ( emmitter_id )
  for(BranchingList::const_iterator cit = branchings.lower_bound( abs(emmitter_id) );
      cit != branchings.upper_bound( abs(emmitter_id) ); ++cit ) {
    IdList ids = cit->second.second;
    if( abs( ids[0] ) == abs( emmitter_id ) ) {
      if( abs(ids[1]) == abs(part_i->id()) && 
	  abs(ids[2]) == abs(part_j->id()) ) {
	return cit->second.first;
      }
      if( abs( ids[1] ) == abs( part_j->id() ) && 
	  abs( ids[2] ) == abs( part_i->id() ) ) {
	swap( part_i, part_j );
	return cit->second.first;
      }
    }
  }
  return SudakovPtr();
}
void PowhegHandler::fillProtoTrees( map< ShowerParticlePtr, HardBranchingPtr > branchingMap, 
				    ProtoTreePtr currentProtoTree ){
  // SHOULD CHECK WITH ME class BUT MOST CONCIVABLE ONES ARE 2->2
  if(branchingMap.size()==4) return;
  for( map< ShowerParticlePtr, HardBranchingPtr >::iterator ita = branchingMap.begin();
       ita != branchingMap.end(); ++ita ) {
    for( map< ShowerParticlePtr, HardBranchingPtr >::iterator itb = branchingMap.begin();
	 itb != ita; ++itb) {
      if( ita->second->status() == HardBranching::Incoming && 
	  itb->second->status() == HardBranching::Incoming ) continue;
      bool incoming = 
	ita->second->status() == HardBranching::Incoming ||
	itb->second->status() == HardBranching::Incoming;
      // make sure particles are in the right order with incoming first
      ShowerParticlePtr part1 = ita->first;
      ShowerParticlePtr part2 = itb->first;  
      if( itb->second->status() == HardBranching::Incoming ) swap( part1, part2 );
      // Get the HardBranching for the pair - this is done such that only a single
      // clustering of each pair of hardBranchings is created.
      HardBranchingPtr currentBranching = getCluster( make_pair( part1, part2 ), 
						      branchingMap, incoming );
      if( ! currentBranching ) continue;
      // branching allowed so make a new Tree out of these branchings
      set< HardBranchingPtr > newTreeBranchings = currentProtoTree->getBranchings();
      ProtoTreePtr newProtoTree = new_ptr( ProtoTree( newTreeBranchings ) );
      bool test1 = newProtoTree->removeBranching( ita->second );
      bool test2 = newProtoTree->removeBranching( itb->second );
      assert(test1 && test2);
      newProtoTree->addBranching( currentBranching );
      map< ShowerParticlePtr, HardBranchingPtr > newBranchingMap = branchingMap;
      assert( newBranchingMap.find( part1 ) != newBranchingMap.end() );
      newBranchingMap.erase( part1 );
      assert( newBranchingMap.find( part2 ) != newBranchingMap.end() );
      newBranchingMap.erase( part2 );
      newBranchingMap.insert( make_pair( currentBranching->branchingParticle(),
					 currentBranching ) );
    
      if( ! repeatProtoTree( newProtoTree ) ) _proto_trees.insert( newProtoTree );
     
      // remove the current tree if it hasn't already been removed
      if( _proto_trees.find( currentProtoTree ) != _proto_trees.end() )
	_proto_trees.erase( currentProtoTree );
      //do recursion
      fillProtoTrees( newBranchingMap, newProtoTree );
    }
  }
}

HardBranchingPtr PowhegHandler::getCluster( pair< ShowerParticlePtr, ShowerParticlePtr > clusterPair,
					    map< ShowerParticlePtr, HardBranchingPtr > theParticles,
					    bool incoming ){
  //look for the clustered pair in _all_clusters
  for( map< HardBranchingPtr , pair< ShowerParticlePtr, ShowerParticlePtr > >::const_iterator 
	 cit = _all_clusters.begin(); cit != _all_clusters.end(); ++cit ){
    if( ( cit->second.first == clusterPair.first && cit->second.second == clusterPair.second ) ||
	( cit->second.first == clusterPair.second && cit->second.second == clusterPair.first ) ){
      return cit->first;
    }
  }
  //what happens if the branching is a cc branching??
  BranchingElement theBranching;
  if( !incoming ) theBranching = allowedFinalStateBranching( clusterPair );
  else theBranching = allowedInitialStateBranching( clusterPair );

  //if branching is not allowed return null hardbranching
  if( !theBranching.first ) return HardBranchingPtr();

  tcPDPtr particle_data;
  if( !incoming ) particle_data = getParticleData( theBranching.second[0] );
  else particle_data = getParticleData( theBranching.second[1] );

  //create clustered hardBranching
  HardBranchingPtr clusteredBranch;
  if( !incoming ){
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass( 0. * MeV );
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );
    clusteredBranch = new_ptr( HardBranching( clustered, theBranching.first,
					      HardBranchingPtr(), HardBranching::Outgoing ) );
  }
  else{
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() - 
      clusterPair.second->momentum();
    pairMomentum.setMass( 0. * MeV );
    //clean up this logic
    if( particle_data->CC() &&
	( clusterPair.first ->id() != theBranching.second[0] ||
	  clusterPair.second->id() != theBranching.second[2] ) ) {
      particle_data = particle_data->CC();
    }
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, false ) );
    clustered->set5Momentum( pairMomentum );
    clusteredBranch = new_ptr( HardBranching( clustered, SudakovPtr(),
					      HardBranchingPtr(), HardBranching::Incoming ) );
    //add back sudakovPtr
    clusteredBranch->backSudakov( theBranching.first );
  }
  _all_clusters.insert( make_pair( clusteredBranch, clusterPair ) );
  //set children relations 
   if( !incoming ){
     clusteredBranch->addChild( theParticles.find( clusterPair.first )->second );	    
     clusteredBranch->addChild( theParticles.find( clusterPair.second )->second );  
   }
   else{
     clusteredBranch->addBackChild( theParticles.find( clusterPair.first )->second );	    
     clusteredBranch->addBackChild( theParticles.find( clusterPair.second )->second );  
   }
  return clusteredBranch;
}

bool PowhegHandler::repeatProtoTree( ProtoTreePtr currentProtoTree ) {
  // loop over all prototrees and see 
  // how many hardbranchings of curentProtoTree are found in each
  for( set< ProtoTreePtr >::const_iterator cit = _proto_trees.begin();
       cit != _proto_trees.end(); ++cit ) {
    unsigned int no_matches = 0;
    for( set< HardBranchingPtr >::const_iterator ckt 
	   = currentProtoTree->getBranchings().begin(); 
	 ckt != currentProtoTree->getBranchings().end(); ckt++ ) {
      if( (*cit)->getBranchings().find( *ckt ) != (*cit)->getBranchings().end() )
	++no_matches;
    }
    // DOES THIS EVER RETURN TRUE ????
    if( no_matches == currentProtoTree->getBranchings().size() ) {
      //      assert(false);
      return true;
    }
  }
  return false;
}

PotentialTree PowhegHandler::doClustering() { 
  assert( _matrixElement );
  // get particles from the XComb object
  ParticleVector outgoing  = lastXCombPtr()->subProcess()->outgoing();
  PPair incoming = lastXCombPtr()->subProcess()->incoming();
  // clear storage 
  _all_clusters.clear();
  _proto_trees.clear();
  _noShowerHists = false;
  // remove backChildRelations from all hardtrees 
  // - to avoid cyclic pointer relations causing memory leaks
  for( size_t ix = 0; ix < _rejectedCKKWTrees.size(); ++ix )
    _rejectedCKKWTrees[ix].first.tree->removeBackChildren();
  for( size_t ix = 0; ix < _hardTrees.size(); ix++ )
    _hardTrees[ix].first.tree->removeBackChildren();
  // clear the hard trees
  _hardTrees.clear();
  _rejectedCKKWTrees.clear();
  map <ShowerParticlePtr,HardBranchingPtr> branchingMap;
  // loops through the FS particles and create hardBranchings
  for( unsigned int i = 0; i < outgoing.size(); ++i) {
    ShowerParticlePtr currentParticle = new_ptr( ShowerParticle( *outgoing[i], 1, true, false ) );
    HardBranchingPtr currentBranching = new_ptr( HardBranching( currentParticle, SudakovPtr(),
								HardBranchingPtr(), 
								HardBranching::Outgoing ) ); 
    branchingMap.insert( make_pair( currentParticle, currentBranching ) );
  }
  // add IS hardBranchings
  ShowerParticlePtr currentParticle = new_ptr( ShowerParticle( *incoming.first, 1, false, false ) );
  HardBranchingPtr currentBranching =  new_ptr( HardBranching( currentParticle, SudakovPtr(),
							       HardBranchingPtr(),
							       HardBranching::Incoming) );
  branchingMap.insert( make_pair( currentParticle, currentBranching ) );  
  currentParticle = new_ptr( ShowerParticle( *incoming.second, 1, false, false ) );
  currentBranching =  new_ptr( HardBranching( currentParticle, SudakovPtr(), HardBranchingPtr(),
					      HardBranching::Incoming ) );
  branchingMap.insert( make_pair( currentParticle, currentBranching ) );

  //create and initialise the first tree
  ProtoTreePtr initialProtoTree = new_ptr( ProtoTree() );
  for( map<ShowerParticlePtr, HardBranchingPtr>::iterator ita = branchingMap.begin();
       ita != branchingMap.end(); ita++ ){
    initialProtoTree->addBranching( ita->second );
  }
  _proto_trees.insert( initialProtoTree );
  //fill _proto_trees with all possible trees
  fillProtoTrees( branchingMap, initialProtoTree );

  double totalWeight = 0.;
  // create a hardtree from each proto tree and fill _hardTrees with angular ordered configs
  for( set< ProtoTreePtr >::const_iterator cit = _proto_trees.begin(); 
       cit != _proto_trees.end(); ++cit ){
    vector<HardBranchingPtr> theBranchings, spaceBranchings;
    (**cit).fixFwdBranchings();
    for( set< HardBranchingPtr >::const_iterator cjt = (*cit)->getBranchings().begin(); 
	 cjt != (*cit)->getBranchings().end(); ++cjt ) {
      theBranchings.push_back( *cjt );
      if( (*cjt)->status() == HardBranching::Incoming ) {
	tHardBranchingPtr spaceLike = *cjt;
	while( spaceLike->parent() ) {
	  spaceLike = spaceLike->parent();
	}
	spaceBranchings.push_back( spaceLike );
      }
    }
    // Create the hard tree
    PotentialTree newTree;
    newTree.tree = new_ptr( CKKWTree( theBranchings, spaceBranchings ,
				      ShowerInteraction::QCD) );
    // check the created CKKWTree corresponds to an allowed LO configuration
    // (does matrix element have a corresponding diagram)
    getDiagram( newTree );
    if( !newTree.diagram ) continue;
    // set the beam particles
    setBeams( newTree.tree );
    //do inverse momentum reconstruction
    if( !evolver()->showerModel()->kinematicsReconstructor()
    	->deconstructHardJets( newTree.tree, evolver(), ShowerInteraction::QCD ) ) {
      generator()->log()<<"\n\nproblem with deconstructhardjets\n\n\n";
      continue;
    }
    newTree.tree->fixFwdBranchings();
    newTree.tree->findNodes();
    if( newTree.tree->checkHardOrdering() && 
	newTree.tree->checkXOrdering() ) {
      //find the wgt and fill _hardTrees map
      double treeWeight = sudakovWeight( newTree.tree );
      _hardTrees.push_back( make_pair( newTree, treeWeight ) );
      totalWeight += treeWeight;
    }
    else {
      _rejectedCKKWTrees.push_back( make_pair( newTree, 1. ) );
      totalWeight += 1.;
    }
  }
  if( _hardTrees.empty() ){
    _unorderedEvents++;
    if( _rejectNonAO  ) return PotentialTree();
  }
  PotentialTree chosen_hardTree;  
  //choose a hardTree from shower probability
  if( _clusterOption == 0 ){
    long treeIndex;
    do{
      treeIndex = UseRandom::irnd( _hardTrees.size() ); 
    } 
    while ( _hardTrees[ treeIndex ].second / totalWeight < UseRandom::rnd() );
    chosen_hardTree.tree = _hardTrees[ treeIndex ].first.tree;
  }
  //choose hardtree with lowest pt
  else if( _clusterOption == 2 ){
    Energy min_pt = Constants::MaxEnergy;
    if( !_hardTrees.empty() ){
      for(unsigned int ix = 0; ix < _hardTrees.size(); ix++ ){
	if( _hardTrees[ix].first.tree->totalPt() < min_pt ) {
	  min_pt = _hardTrees[ix].first.tree->totalPt();
	  chosen_hardTree.tree = _hardTrees[ix].first.tree;
	}
      }
    }
    else{
      for(unsigned int ix = 0; ix < _rejectedCKKWTrees.size(); ix++ ){
	if( _rejectedCKKWTrees[ix].first.tree->totalPt() < min_pt ){
	  min_pt = _rejectedCKKWTrees[ix].first.tree->totalPt();
	  chosen_hardTree.tree = _rejectedCKKWTrees[ix].first.tree;
	}
      }
    }
  }
  //choose hardtree with highest prob (cluster option = 3)
  else if( _clusterOption == 3 ) {
    double max_prob = 0.;
    for(unsigned int ix = 0; ix < _hardTrees.size(); ix++ ){
      if( _hardTrees[ ix ].second > max_prob ){
	max_prob = _hardTrees[ ix ].second;
	chosen_hardTree.tree = _hardTrees[ix].first.tree;
      }
      if( _hardTrees[ ix ].second == max_prob && UseRandom::rndbool() ){
	max_prob = _hardTrees[ ix ].second;
	chosen_hardTree.tree = _hardTrees[ix].first.tree;
      }
    }
  }
  else {
    assert(false);
  }
  // re-do momentum deconstruction (has been overridden by other trees otherwise)
  if(! chosen_hardTree.tree ) {
    if( ! _rejectNoHistories ) _noShowerHists = true;
    return PotentialTree();
  }
  assert( chosen_hardTree.tree );
  chosen_hardTree.tree->fixFwdBranchings();
  setBeams( chosen_hardTree.tree );
  getDiagram( chosen_hardTree);
  bool result = evolver()->showerModel()->kinematicsReconstructor()
    ->deconstructHardJets( chosen_hardTree.tree, evolver(),ShowerInteraction::QCD );
  assert( result );
  assert( chosen_hardTree.tree );
  ++_trees_created;
  return chosen_hardTree;
}

void PowhegHandler::fixColours(tPPtr parent, tPPtr child1, tPPtr child2) {
  // the different possible cases
  if(parent->dataPtr()->iColour()==PDT::Colour3&&
     child1->dataPtr()->iColour()==PDT::Colour3&&
     child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->colourLine()->addColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    temp->addColoured(child1);
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3&&
	  child2->dataPtr()->iColour()==PDT::Colour3&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->colourLine()->addColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    temp->addColoured(child2);
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    temp->addColoured(child2);
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    temp->addColoured(child1);
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour8&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    if(UseRandom::rndbool(0.5)) {
      child1->colourLine()->addColoured(parent);
      child2->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child1->antiColourLine();
      temp->addColoured(child2);
      child2->colourLine()->join(temp);
    }
    else {
      child2->colourLine()->addColoured(parent);
      child1->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child2->antiColourLine();
      temp->addColoured(child1);
      child1->colourLine()->join(temp);
    }
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour3&&
	  child2->dataPtr()->iColour()==PDT::Colour3bar) {
    child1->colourLine()->addColoured(parent);
    child2->antiColourLine()->addAntiColoured(parent);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour3) {
    child2->colourLine()->addColoured(parent);
    child1->antiColourLine()->addAntiColoured(parent);
  }  
  else {
    assert(false);
  }
}

BranchingElement PowhegHandler::
allowedFinalStateBranching( pair< ShowerParticlePtr, ShowerParticlePtr > & clusterPair ) {
  // check with normal ID's
  pair< long, long > ptest = make_pair( clusterPair.first->id(), clusterPair.second->id() );
  map< pair< long, long >, pair< SudakovPtr, IdList > >::const_iterator 
    split = _allowedFinal.find(ptest);
  if( split != _allowedFinal.end() ) {
    if(  split->second.second[1] != ptest.first ) swap( clusterPair.first, clusterPair.second );
    return split->second;
  }
  // check with CC
  if( clusterPair.first->dataPtr()->CC() ) ptest.first  *= -1;
  if( clusterPair.second->dataPtr()->CC() ) ptest.second *= -1;
  split = _allowedFinal.find( ptest );
  if( split != _allowedFinal.end() ) {
    //cc the idlist
    //this will only be for qbar g clusterings
    BranchingElement ccBranch = split->second;
    if( getParticleData( ccBranch.second[0] )->CC() ) ccBranch.second[0] *= -1;
    if( getParticleData( ccBranch.second[1] )->CC() ) ccBranch.second[1] *= -1;
    if( getParticleData( ccBranch.second[2] )->CC() ) ccBranch.second[2] *= -1;
    if( split->second.second[1] !=  ptest.first ) swap( clusterPair.first, clusterPair.second );
    return ccBranch;
  }
  // not found found null pointer
  return make_pair( SudakovPtr(), IdList() );
}

BranchingElement PowhegHandler::
allowedInitialStateBranching( pair< ShowerParticlePtr, ShowerParticlePtr > & br ) {
  // veto top
  if( abs( br.first->id() ) == ParticleID::t ||
      abs( br.second->id() ) == ParticleID::t )
    return make_pair( SudakovPtr(), IdList() );
  //is initial parton an antiparticle
  bool cc = br.first->id() < 0;
  //gives range of _allowedInitial with matching first abs( id )
  pair< multimap< long, pair< SudakovPtr, IdList > >::const_iterator,
    multimap< long, pair< SudakovPtr, IdList > >::const_iterator >
    location = _allowedInitial.equal_range( abs( br.first->id() ) );
  //iterates over this range
  for( multimap< long, pair< SudakovPtr, IdList> >::const_iterator it = location.first;
       it != location.second; ++it ) {
    //test id for second particle in pair
    long idtest = it->second.second[2];
    //if it is antiparticle *= -1
    if( cc && getParticleData( idtest )->CC() ) idtest *= -1;
    //does second id match the test
    if( idtest == br.second->id() ) return it->second;
    //if the the IS parton is a gluon and charge conjugate of second parton mathes accept
    if( idtest == -br.second->id() &&
        ! br.first->dataPtr()->CC() ) return it->second;
  }
  // not found found null pointer
  return make_pair(SudakovPtr(),IdList());
}

namespace {
struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

void PowhegHandler::getDiagram(PotentialTree & tree) {
  if(tree.diagrams.empty()) {
    set<HardBranchingPtr>::const_iterator cit;
    tcPDPair incoming;
    multiset<tcPDPtr,ParticleOrdering> outgoing;  
    //get the incoming and outgoing partons involved in hard process
    for( cit = tree.tree->branchings().begin(); 
	 cit != tree.tree->branchings().end(); ++cit ){ 
      if( (*cit)->status() ==HardBranching::Incoming) {
	HardBranchingPtr parent = *cit;
	while(parent->parent()) parent = parent->parent();
	if( parent->branchingParticle()->momentum().z()>ZERO )
	  incoming.first  = (*cit)->branchingParticle()->dataPtr();
	else
	  incoming.second = (*cit)->branchingParticle()->dataPtr();
      }
      else {
	outgoing.insert( (*cit)->branchingParticle()->dataPtr() );
      }
    }
    if(!incoming.first || !incoming.second) return;
    pair<string,string> tag;
    tag.first  = incoming.first  ->PDGName() + "," + incoming.second->PDGName() + "->";
    tag.second = incoming.second ->PDGName() + "," + incoming.first ->PDGName() + "->";

    string tag_out;
    for ( multiset<tcPDPtr,ParticleOrdering>::iterator i = outgoing.begin();
	  i != outgoing.end(); ++i ) {
      if ( i != outgoing.begin() ) tag_out += ",";
      tag_out += (**i).PDGName();
    }
    tag.first  += tag_out;
    tag.second += tag_out;
    // find the diagrams
    for ( int i = 0, N = _matrixElement->diagrams().size(); i < N; ++i ) {
      string temp = _matrixElement->diagrams()[i]->getTag();
      if ( temp == tag.first || temp == tag.second )
	tree.diagrams.push_back(_matrixElement->diagrams()[i]);
    }
  }
  if(tree.diagrams.empty()) return;
  // construct a set of on-shell momenta for the hard collison
  vector<Lorentz5Momentum> meMomenta;
  vector<tcPDPtr> mePartonData;
  PVector particles;
  set<HardBranchingPtr>::const_iterator it;
  // incoming particles
  for( it = tree.tree->branchings().begin();it != tree.tree->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Incoming ) {
      meMomenta.push_back( (**it).branchingParticle()->momentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
      particles.push_back( (**it).branchingParticle() );
    }
  }
  assert(particles.size()==2);
  for( it = tree.tree->branchings().begin(); it != tree.tree->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Outgoing ) {
      meMomenta.push_back( (**it).branchingParticle()->momentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
      particles.push_back( (**it).branchingParticle() );
    }
  }
  const cPDVector partons = tree.diagrams[0]->partons();
  // order of the incoming partons
  if(mePartonData[0] != partons[0]) {
    swap( mePartonData[0], mePartonData[1] );
    swap( meMomenta[0], meMomenta[1] );
    swap( particles[0], particles[1] );
  }
  // order of the outgoing partons
  for(unsigned int ix=2;ix<partons.size();++ix) {
    for(unsigned int iy=ix;iy<meMomenta.size();++iy) {
      if(partons[ix]==mePartonData[iy]) {
	if(ix!=iy) {
	  swap(mePartonData[ix],mePartonData[iy]);
	  swap(meMomenta[ix],meMomenta[iy]);
	  swap(particles[ix],particles[iy]);
	}
	break;
      }
    }
  }
  // boost to the CMF frame
  Lorentz5Momentum prest(meMomenta[0]+meMomenta[1]);
  LorentzRotation R(-prest.boostVector());
  // and then to put beams along the axis
  Lorentz5Momentum ptest = R*meMomenta[0];
  Axis axis(ptest.vect().unit());
  if(axis.perp2()>0.) {
    R.rotateZ(-axis.phi());
    R.rotateY(-acos(axis.z()));
  }
  for( unsigned int ix = 0; ix < meMomenta.size(); ++ix )
    meMomenta[ix].transform(R);
  // now rescale to put on shell
  Energy Ebeam = 0.5 * ( meMomenta[0].e() + meMomenta[1].e() );
  for( unsigned int i = 0; i < 2; ++i ) {
    meMomenta[i].setZ( meMomenta[i].z() / abs(meMomenta[i].z()) * Ebeam  );
    meMomenta[i].setT( Ebeam );
  }
  Energy2 s = 4.0 * sqr(Ebeam);
  Energy m1 = mePartonData[2]->mass();
  Energy m2 = mePartonData[3]->mass();
  double lambda = 0.25/Ebeam/meMomenta[2].rho() * 
    sqrt( ( s - sqr(m1+m2) ) * ( s - sqr(m1-m2) ) );
  for( unsigned int i = 2; i < meMomenta.size(); ++i ) {
    meMomenta[i] *= lambda;
    meMomenta[i].setMass(mePartonData[i]->mass());
    meMomenta[i].rescaleEnergy();
  }
  // incoming pair
  PPair in( mePartonData[0]->produceParticle( meMomenta[0] ),
	    mePartonData[1]->produceParticle( meMomenta[1] ) );
  // outgoing
  PVector out;
  for(unsigned int ix=2;ix<meMomenta.size();++ix) {
    out.push_back(mePartonData[ix]->produceParticle(meMomenta[ix]));
  }
  // call the matrix element to initialize
  _matrixElement->setKinematics(in,out);
  _matrixElement->dSigHatDR();
  // select the diagram
  if(!tree.diagram)
    tree.diagram = tree.diagrams[_matrixElement->diagram(tree.diagrams)];
  // get the colour structure
  if(!tree.cl) {
    Selector<const ColourLines *> sel = _matrixElement->colourGeometries(tree.diagram);
    tree.cl = sel.select(rnd());
  }
  PVector slike;
  tPVector ret;
  slike.push_back(particles[0]);
  Ptr<Tree2toNDiagram>::pointer diagram2 = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::pointer>(tree.diagram);
  for ( int i = 1; i < diagram2->nSpace() - 1; ++i )
    slike.push_back(diagram2->allPartons()[i]->produceParticle());
  slike.push_back(particles[1]);
  ret = tPVector(slike.begin(), slike.end());
  int io = particles.size();
  PVector tlike(diagram2->allPartons().size() - diagram2->nSpace());
  for ( int i = diagram2->allPartons().size() - 1; i >=  diagram2->nSpace(); --i ) {
    int it = i - diagram2->nSpace();
    pair<int,int> ch = diagram2->children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = particles[--io];
    } 
    else {
      Lorentz5Momentum p = tlike[ch.first - diagram2->nSpace()]->momentum() +
 	tlike[ch.second - diagram2->nSpace()]->momentum();
      tlike[it] = diagram2->allPartons()[i]->produceParticle(p);
    }
  }
  ret.insert( ret.end(), tlike.begin(), tlike.end() );
  tree.cl->connect(ret);
  for( unsigned int ix = 0; ix < ret.size(); ++ix ) {
    PVector::iterator it = find( particles.begin(), particles.end(),ret[ix] );
    if( it == particles.end() ) {
      ColinePtr line = ret[ix]->colourLine();
      if(line) line->removeColoured(ret[ix]);
      line = ret[ix]->antiColourLine();
      if(line) line->removeAntiColoured(ret[ix]);
    }
  }
  // now the colours of the rest of the particles
  for( set<HardBranchingPtr>::const_iterator it = tree.tree->branchings().begin();
       it!=tree.tree->branchings().end(); ++it ) (**it).fixColours();
  // now make the colour connections in the tree
  ShowerParticleVector branchingParticles;
  map<ShowerParticlePtr,HardBranchingPtr> branchingMap;
  for( set< HardBranchingPtr >::iterator it = tree.tree->branchings().begin();
       it != tree.tree->branchings().end(); ++it ) {
    branchingParticles.push_back((**it).branchingParticle());
    branchingMap.insert(make_pair((**it).branchingParticle(),*it));
  }
  // find the colour partners
  evolver()->showerModel()->partnerFinder()
    ->setInitialEvolutionScales(branchingParticles,false,
				ShowerInteraction::QCD,true);
  for(unsigned int ix=0;ix<branchingParticles.size();++ix) {
    if(branchingParticles[ix]->partner()) {
      HardBranchingPtr partner = branchingMap[branchingParticles[ix]->partner()];
      branchingMap[branchingParticles[ix]]->colourPartner(partner);
    }
  }
}

double PowhegHandler::splittingFnWeight( CKKWTreePtr theTree ){
  double splitFnWgt = 1.;
  for( map< HardBranchingPtr, Energy>::const_iterator cit = theTree->getNodes().begin();
       cit != theTree->getNodes().end(); ++cit ) {
    vector< long > ids;
    ids.push_back( cit->first->branchingParticle()->id() );
    assert( ! cit->first->children().empty() );
    ids.push_back( cit->first->children()[0]->branchingParticle()->id() );
    ids.push_back( cit->first->children()[1]->branchingParticle()->id() );
    double z = cit->first->children()[0]->z();
    Energy2 t = z * ( 1. - z ) * sqr( cit->second );
    splitFnWgt *= cit->first->sudakov()->splittingFn()->P( z, t, ids, true );
  }
  return splitFnWgt;
}

double PowhegHandler::sudakovWeight( CKKWTreePtr theTree ) {
  //ktcut for sudakovs
  Energy kt_cut = _highestMult ? 
    theTree->lowestPt( _jetMeasureMode, _s ) : _showerVeto->getVeto();
  //alphaS weight
  double alphaWgt = 1.;
  for( map< HardBranchingPtr, Energy >::const_iterator cit = theTree->getNodes().begin(); 
       cit != theTree->getNodes().end(); ++cit ) {
    if( cit->first->status() == HardBranching::Incoming ) {
      double wgt =  _alphaS->value( sqr( cit->first->pT() ) ) / _alphaSMG;
      assert( !isnan(wgt) );
      alphaWgt *= wgt;
    }
    else if( ! cit->first->children().empty() ) {
      double wgt = _alphaS->value( sqr( cit->first->children()[0]->pT() ) ) / _alphaSMG;
      assert( !isnan(wgt) );
      alphaWgt *= wgt;
    }
    else {
      assert(false);
    }
  }
  assert(!isnan(alphaWgt) );
  double PDFweight = pdfWeight( theTree );
  assert( !isnan( PDFweight ) );
  return alphaWgt*PDFweight;
}

void PowhegHandler::setBeams( CKKWTreePtr tree ) {
  PPair beams = lastXCombPtr()->lastParticles();
  //remove children of beams
  PVector beam_children = beams.first->children();
  for( unsigned int ix = 0; ix != beam_children.size(); ix++ )
    beams.first->abandonChild( beam_children[ix] );
  beam_children = beams.second->children();
  for( unsigned int ix = 0; ix != beam_children.size(); ix++ )
    beams.second->abandonChild( beam_children[ix] );
  if( (**tree->incoming().begin()).branchingParticle()->momentum().z() /
      beams.first->momentum().z() < 0.)
    swap( beams.first, beams.second );
  set<HardBranchingPtr>::iterator it = tree->incoming().begin();
  HardBranchingPtr br = *it;
  br->beam( beams.first );
  while ( !br->children().empty() ) {
    for(unsigned int ix = 0; ix < br->children().size(); ++ix ) {
      if( br->children()[ix]->status() == HardBranching::Incoming ) {
	br = br->children()[ix];
	break;
      }
    }
    br->beam( beams.first );
  }
  ++it;
  br = *it;
  br->beam( beams.second );
  while ( !br->children().empty() ) {
    for( unsigned int ix = 0; ix < br->children().size(); ++ix ) {
      if( br->children()[ix]->status() == HardBranching::Incoming ) {
	br = br->children()[ix];
	break;
      }
    }
    br->beam( beams.second );
  }
}

bool PowhegHandler::updateSubProcess() {
  assert( hardTree_.diagram && hardTree_.tree );
  Ptr<Tree2toNDiagram>::pointer diagram = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::pointer>(hardTree_.diagram);
  assert(diagram);
  set<HardBranchingPtr>::const_iterator it; 
  map< ColinePtr, ColinePtr> colourMap;
  PPair incoming;
  // loop over the branchings and sort out incoming particles
  for( it = hardTree_.tree->branchings().begin();
       it != hardTree_.tree->branchings().end(); ++it) {
    if( (*it)->status() == HardBranching::Outgoing ) continue;
    PPtr newParticle = new_ptr( Particle( (**it).branchingParticle()->dataPtr() ) );
    newParticle->set5Momentum( (**it).showerMomentum() );
    if( (**it).branchingParticle()->colourLine() ) {
      map< ColinePtr, ColinePtr>::iterator loc 
	= colourMap.find( (**it).branchingParticle()->colourLine() );
      if( loc != colourMap.end() ) {
	loc->second->addColoured( newParticle );
      }
      else {
	ColinePtr newLine = new_ptr( ColourLine() );
	colourMap[ (**it).branchingParticle()->colourLine() ] = newLine;
	newLine->addColoured( newParticle );
      }
    }
    if( (**it).branchingParticle()->antiColourLine() ) {
      map< ColinePtr, ColinePtr> ::iterator loc 
	= colourMap.find( (**it).branchingParticle()->antiColourLine() );
      if( loc != colourMap.end() ) {
	loc->second->addAntiColoured( newParticle );
      }
      else {
	ColinePtr newLine = new_ptr( ColourLine() );
	colourMap[ (**it).branchingParticle()->antiColourLine() ] = newLine;
	newLine->addAntiColoured( newParticle );
      }
    }
    if( lastXCombPtr()->subProcess()->incoming().first->momentum().z() /
	newParticle->momentum().z() > 0. )
      incoming.first  = newParticle;
    else
      incoming.second = newParticle;
  }
  bool mirror = 
    incoming.first ->dataPtr() != diagram->partons()[0] &&
    incoming.second->dataPtr() != diagram->partons()[1];
  // create the new subprocess
  SubProPtr newSubProcess =
    new_ptr( SubProcess( incoming, lastXCombPtr()->subProcess()->collision(),
			 lastXCombPtr()->subProcess()->handler() ) );
  // add the spacelike intermediates
  PVector slike;
  slike.push_back( !mirror ? incoming.first : incoming.second);
  for ( int i = 1; i < diagram->nSpace() - 1; ++i )
    slike.push_back(diagram->allPartons()[i]->produceParticle());
  slike.push_back( !mirror ? incoming.second : incoming.first);
  tPVector ret = tPVector(slike.begin(), slike.end());
  for ( size_t i = 1; i < slike.size() - 1; ++i ) {
    slike[i-1]->addChild(slike[i]);
    newSubProcess->addIntermediate(slike[ mirror ? i: slike.size() - 1 - i], false);
  }
  // get the parton bins from the parton extractor if first time
  static bool first = true;
  if(first) {
    first = false;
    Energy e1 = lastXCombPtr()->lastParticles().first ->momentum().e();
    Energy e2 = lastXCombPtr()->lastParticles().second->momentum().e();
    Energy emax = 2.0*sqrt(e1*e2);
    cPDPair inData = make_pair(lastXCombPtr()->lastParticles().first ->dataPtr(),
			       lastXCombPtr()->lastParticles().second->dataPtr());
    cuts_->initialize(sqr(emax),0.5*log(e1/e2));
    partonBins_ = partonExtractor_->getPartons(emax, inData, *cuts_);
  }
  // get the parton bins for this event
  tcPBPair sel;
  for ( int i = 0, N = partonBins_.size(); i < N; ++i ) {
    tcPBPtr bin = partonBins_[i].first;
    tPPtr p = incoming.first;
    while ( bin && p ) {
      if ( p->dataPtr() != bin->parton() ) break;
      bin = bin->incoming();
      p = p != lastXCombPtr()->lastParticles().first ?
	lastXCombPtr()->lastParticles().first : PPtr();
    }
    if ( bin || p ) continue;
    bin = partonBins_[i].second;
    p = incoming.second;
    while ( bin && p ) {
      if ( p->dataPtr() != bin->parton() ) break;
      bin = bin->incoming();
      p = p != lastXCombPtr()->lastParticles().second ?
	lastXCombPtr()->lastParticles().second : PPtr();
    }
    if ( bin || p ) continue;
    sel = partonBins_[i];
    break;
  }
  if ( !sel.first || !sel.second ) Throw<Exception>()
				     << "Could not find appropriate "
				     << "PartonBin objects for event in "
				     << "PowhegHandler " << Exception::runerror;
  // create the new parton bin instances
  Direction<0> dir(true);
  PBIPair partonBinInstances;
  // temporary mother/daugther settings
  // to get parton bin instances correct
  lastXCombPtr()->lastParticles().first ->addChild(incoming.first );
  lastXCombPtr()->lastParticles().second->addChild(incoming.second);
  // make the parton bin instances
  partonBinInstances.first =
    new_ptr(PartonBinInstance(incoming.first, sel.first,
			      lastXCombPtr()->partonBinInstances().first->scale()));
  dir.reverse();
  partonBinInstances.second =
    new_ptr(PartonBinInstance(incoming.second, sel.second,
			      lastXCombPtr()->partonBinInstances().second->scale()));
  // remove temporary mother/daugther settings
  lastXCombPtr()->lastParticles().first ->abandonChild(incoming.first );
  lastXCombPtr()->lastParticles().second->abandonChild(incoming.second);
  // set the parton bin instances
  lastXCombPtr()->setPartonBinInstances(partonBinInstances,
					lastXCombPtr()->lastScale());
  // momenta of the time like partons
  vector<Lorentz5Momentum> meMomenta;
  vector<tcPDPtr> mePartonData;
  for( it = hardTree_.tree->branchings().begin(); 
       it != hardTree_.tree->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Outgoing ) {
      meMomenta.push_back( (**it).showerMomentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
    }
  }
  // order of the outgoing partons
  for(int ix=2;ix<int(diagram->partons().size());++ix) {
    for(int iy=ix-2;iy<int(meMomenta.size());++iy) {
      if(diagram->partons()[ix]==mePartonData[iy]) {
	if(ix!=iy) {
	  swap(mePartonData[ix-2],mePartonData[iy]);
	  swap(meMomenta   [ix-2],meMomenta   [iy]);
	}
	break;
      }
    }
  }
  // time like particles
  int io = meMomenta.size();
  PVector tlike(diagram->allPartons().size() - diagram->nSpace());
  tPVector out;
  for ( int i = diagram->allPartons().size() - 1; i >=  diagram->nSpace(); --i ) {
    int it = i - diagram->nSpace();
    pair<int,int> ch = diagram->children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = diagram->allPartons()[i]->produceParticle(meMomenta[--io]);
    } 
    else {
      Lorentz5Momentum p = tlike[ch.first - diagram->nSpace()]->momentum() +
	tlike[ch.second - diagram->nSpace()]->momentum();
      tlike[it] = diagram->allPartons()[i]->produceParticle(p);
    }
    if ( diagram->parent(i) < diagram->nSpace() ) {
      slike[diagram->parent(i)]->addChild(tlike[it]);
      if ( diagram->parent(i) == diagram->nSpace() - 2 )
	slike[diagram->parent(i) + 1]->addChild(tlike[it]);
    }
    if ( !iso ) {
      tlike[it]->addChild(tlike[ch.first - diagram->nSpace()]);
      tlike[it]->addChild(tlike[ch.second - diagram->nSpace()]);
    }
    if ( iso )  out.push_back(tlike[it]);
    else        newSubProcess->addIntermediate(tlike[it], false);
  }
  ret.insert(ret.end(), tlike.begin(), tlike.end());   
  for ( int i = 0, N = out.size(); i < N; ++i )
    newSubProcess->addOutgoing(out[mirror ? i: out.size() - i - 1], false);
  for ( PVector::size_type i = 0; i < slike.size() - 2; ++i ) {
    pair<int,int> ch = diagram->children(i);
    slike[ch.first]->set5Momentum(slike[i]->momentum() -
				  tlike[ch.second - diagram->nSpace()]->momentum());
  }
  const ColourLines & cl = _matrixElement->selectColourGeometry(hardTree_.diagram);
  cl.connect(ret);
  // set the subprocess
  lastXCombPtr()->subProcess( newSubProcess );  
  return true;
} 

bool ProtoTree::fixFwdBranchings(){
  //loop over the spacelike hardbranchings in hard process
  //not sure if this has cleared all previous set relations
  set< HardBranchingPtr >::iterator it;
  for( it = branchings_.begin(); it != branchings_.end(); ++it ){
    HardBranchingPtr current = *it;
    if( current->status() == HardBranching::Outgoing ) continue;
    current->clearChildren();
    current->sudakov( SudakovPtr() );
    if( current->backChildren().empty() ) continue;
    vector< HardBranchingPtr > backChildren = current->backChildren();
    while( ! backChildren.empty() ){
      assert( backChildren.size() == 2 );
      if( !backChildren[0] || !backChildren[1] ) continue;
      //remove any exiting children
      backChildren[0]->clearChildren();
      backChildren[0]->addChild( current );
      backChildren[0]->addChild( backChildren[1] );
      current->parent( backChildren[0] );
      backChildren[1]->parent ( backChildren[0] );
      assert ( current->backSudakov() );
      backChildren[0]->sudakov( current->backSudakov() );
      //continue along incoming line (always the first child)
      current = backChildren[0];
      backChildren = backChildren[0]->backChildren();
    } 
  }
  return true;
}

double PowhegHandler::pdfWeight( CKKWTreePtr tree) {
  double current_pdf;
  double pdf_wgt = 1.;
  double x_frac;
  for( set< HardBranchingPtr >::const_iterator cit = tree->branchings().begin(); 
       cit != tree->branchings().end(); ++cit ) {
    if( (*cit)->status() !=HardBranching::Incoming ) continue;
    //find the pdf for this line
    HardBranchingPtr currentParticle = *cit;
    tcPDPtr beamData = currentParticle->beam()->dataPtr();
    tcPDFPtr pdf = dynamic_ptr_cast< Ptr<BeamParticleData>::const_pointer>(beamData)->pdf();
    if(!pdf) continue;
    // should also implement pdf freeze scale
    Lorentz5Momentum beamMom = currentParticle->beam()->momentum();
    // this assumes we are in c.o.m frame of beams with beam dirn along z axis
    Lorentz5Momentum otherBeam = beamMom;
    otherBeam.setZ( -otherBeam.z() );
    tcPDPtr theBeam = currentParticle->beam()->dataPtr();
    tcPDPtr theParton = currentParticle->branchingParticle()->dataPtr();
    Energy2 scale = lastXCombPtr()->lastSHat();
    x_frac = currentParticle->x_frac();
    //include born pdf factor (that shower would have had)
    current_pdf = pdf->xfx( theBeam, theParton, scale, x_frac ) / x_frac;
    if( current_pdf <= 0. ) return 0.;
    pdf_wgt *= current_pdf;
    // trace line back with a fraction of pdfs for each node
    while( currentParticle->parent() ) {
      scale = sqr( currentParticle->parent()->scale() );
      theParton = currentParticle->branchingParticle()->dataPtr();
      current_pdf = pdf->xfx( theBeam, theParton, 
			      scale, x_frac ) / x_frac;
      // / x_frac;
      if( current_pdf <= 0. ) return 0.;
      pdf_wgt /= current_pdf;
      currentParticle = currentParticle->parent();
      x_frac = currentParticle->x_frac();
      theParton = currentParticle->branchingParticle()->dataPtr();
      current_pdf = pdf->xfx( theBeam, theParton, 
			      scale, x_frac ) / x_frac;
      // / x_frac;
      if( current_pdf <= 0. ) return 0.;
      pdf_wgt *= current_pdf;
    }
    current_pdf = pdf->xfx( theBeam, theParton, 
			    sqr( _pdfScale ), x_frac ) / x_frac;
    if( current_pdf <= 0. ) return 0.;
    pdf_wgt /= current_pdf;
  }
  return pdf_wgt;
}
