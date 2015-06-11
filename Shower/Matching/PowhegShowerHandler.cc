// -*- C++ -*-
//
// PowhegShowerHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegShowerHandler class.
//

#include <config.h>
#include "PowhegShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

// include theses to have complete types
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/PDF/MPIPDF.h"
#include "Herwig++/PDF/MinBiasPDF.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/PDF/HwRemDecayer.h"

#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/HardBranching.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include "Herwig++/MatrixElement/HwMEBase.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/DiagramBase.fh"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"


using namespace Herwig;

namespace {
struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

IBPtr PowhegShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegShowerHandler::fullclone() const {
  return new_ptr(*this);
}

HardTreePtr PowhegShowerHandler::generateCKKW(ShowerTreePtr) const {
  // hard subprocess
  tSubProPtr sub = lastXCombPtr()->subProcess();
  // real emission sub-process
  tSubProPtr real = Factory()->hardTreeSubprocess();
  // born emitter
  emitter_   = Factory()->hardTreeEmitter();
  spectator_ = Factory()->hardTreeSpectator();
  // if no hard emission return
  if ( !(real && emitter_>-1) )
    return HardTreePtr();
  // check emission
  if(sub->outgoing().size()>=real->outgoing().size())
    return HardTreePtr();

  tStdXCombPtr lastXC = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  tStdXCombPtr headXC = lastXC->head();

  if (headXC)
    matrixElement_ = dynamic_ptr_cast<MEPtr>(headXC->matrixElement());
  else if (lastXC)
    matrixElement_ = dynamic_ptr_cast<MEPtr>(lastXC->matrixElement());

  if (lastXC){
    tStdXCombPtr projector= lastXC->lastProjector();
    if (projector){
      matrixElement_ = dynamic_ptr_cast<MEPtr>(projector->matrixElement());
      setSubtractionIntegral(true);
    }
    else 
      setSubtractionIntegral(false);
  }
  assert(matrixElement_);

  // create a hard tree by clustering the event
  hardTree(doClustering(real));
  // Get the HardTree from the CKKW handler.
  CKKWTreePtr hardtree = hardTree().tree();
  // zero to avoid MPI problems
  Factory()->setHardTreeEmitter(-1);
  Factory()->setHardTreeSubprocess(SubProPtr());
  return hardtree;
}

PotentialTree PowhegShowerHandler::doClustering(tSubProPtr real) const {
  // clear storage of the protoTrees
  protoBranchings().clear();
  protoTrees().clear();
  hardTrees_.clear();
  assert( matrixElement() );
  // extract the XComb for the Born process
  tStdXCombPtr lastXC;  
  if (subtractionIntegral()){
    tStdXCombPtr lastXCReal = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
    lastXC = lastXCReal->lastProjector();
  }
  else
    lastXC = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  const StandardXComb xc= *lastXC;
  // get the particles for the born process
  PPair          incomingBorn = xc.subProcess()->incoming();
  ParticleVector outgoingBorn = xc.subProcess()->outgoing();
  // get particles from the XComb object for the  real process
  ParticleVector outgoing  = real->outgoing();
  PPair incoming           = real->incoming();
  // loop through the FS particles and create ProtoBranchings
  for( unsigned int i = 0; i < outgoing.size(); ++i) {
    tPPtr parent = outgoing[i]->parents()[0];
    ProtoBranchingPtr currentBranching =
      new_ptr(ProtoBranching(outgoing[i]->dataPtr(),HardBranching::Outgoing,
			     outgoing[i]->momentum(),tSudakovPtr()));
    currentBranching->    colourLine(outgoing[i]->    colourLine());
    currentBranching->antiColourLine(outgoing[i]->antiColourLine());
    protoBranchings().insert(currentBranching);
  }
  // add IS hardBranchings
  ProtoBranchingPtr currentBranching = 
    new_ptr(ProtoBranching(incoming.first ->dataPtr(),HardBranching::Incoming,
			   incoming.first ->momentum(),tSudakovPtr()));
  currentBranching->    colourLine(incoming.first->    colourLine());
  currentBranching->antiColourLine(incoming.first->antiColourLine());
  protoBranchings().insert(currentBranching);
  currentBranching =
    new_ptr(ProtoBranching(incoming.second->dataPtr(),HardBranching::Incoming,
			   incoming.second->momentum(),tSudakovPtr()));
  currentBranching->    colourLine(incoming.second->    colourLine());
  currentBranching->antiColourLine(incoming.second->antiColourLine());
  protoBranchings().insert(currentBranching);
  // create and initialise the first tree
  ProtoTreePtr initialProtoTree = new_ptr( ProtoTree() );
  for(set<ProtoBranchingPtr>::const_iterator it=protoBranchings().begin();
      it!=protoBranchings().end();++it) { 
    initialProtoTree->addBranching(*it);
  }
  // fill _proto_trees with all possible trees
  protoTrees().insert(initialProtoTree );
  fillProtoTrees( initialProtoTree , xc.mePartonData()[emitter_]->id() );
  // create a HardTree from each ProtoTree and fill hardTrees()
  for( set< ProtoTreePtr >::const_iterator cit = protoTrees().begin(); 
       cit != protoTrees().end(); ++cit ) {
    set<tPPtr> bornParticles(outgoingBorn.begin(),outgoingBorn.end());
    bornParticles.insert(incomingBorn.first );
    bornParticles.insert(incomingBorn.second);
    PotentialTree newTree;
    newTree.tree((**cit).createHardTree());
    // new check based on the colour structure
    map<ColinePtr,ColinePtr> cmap;
    // make the colour connections in the tree
    ShowerParticleVector branchingParticles;
    map<ShowerParticlePtr,HardBranchingPtr> branchingMap;
    bool matched(true);
    int iemitter(-1); 
    HardBranchingPtr emitter;
    map<int,HardBranchingPtr> locMap;
    for( set< HardBranchingPtr >::iterator it = newTree.tree()->branchings().begin();
	 it != newTree.tree()->branchings().end(); ++it ) {
      matched = true;
      // map the particle to the branching for future use
      branchingParticles.push_back((**it).branchingParticle());
      branchingMap.insert(make_pair((**it).branchingParticle(),*it));
      tPPtr bornPartner;
      if((**it).status()==HardBranching::Incoming) {
	HardBranchingPtr parent=*it;
	while(parent->parent()) {
	  parent = parent->parent();
	};
	if(parent->branchingParticle()->momentum().z()/incomingBorn.first->momentum().z()>0.) {
	  bornPartner = incomingBorn.first;
	  if(!parent->children().empty()) {
	    iemitter = 0;
	    emitter = *it;
	  }
	  locMap[0] = *it;
	}
	else {
	  bornPartner = incomingBorn.second;
	  if(!parent->children().empty()) {
	    iemitter = 1;
	    emitter = *it;
	  }
	  locMap[1] = *it;
	}
      }
      else {
	Energy2 dmin( 1e30*GeV2 );
	for(set<tPPtr>::const_iterator bit=bornParticles.begin();bit!=bornParticles.end();
	    ++bit) {
	  if((**it).branchingParticle()->id()!=(**bit).id()) continue;
	  if(*bit==incomingBorn.first||*bit==incomingBorn.second) continue;
	  Energy2 dtest =
	    sqr( (**bit).momentum().x() - (**it).branchingParticle()->momentum().x() ) +
	    sqr( (**bit).momentum().y() - (**it).branchingParticle()->momentum().y() ) +
	    sqr( (**bit).momentum().z() - (**it).branchingParticle()->momentum().z() ) +
	    sqr( (**bit).momentum().t() - (**it).branchingParticle()->momentum().t() );
	  dtest += 1e10*sqr((**bit).momentum().m()-(**it).branchingParticle()->momentum().m());
	  if( dtest < dmin ) {
	    bornPartner = *bit;
	    dmin = dtest;
	  }
	}
	// find the map
	int iloc(-1);
	for(unsigned int ix=0;ix<outgoingBorn.size();++ix) {
	  if(outgoingBorn[ix]==bornPartner) {
	    iloc = ix+2;
	    break;
	  }
	}
	if(!(**it).children().empty()) {
	  emitter = *it;
	  iemitter = iloc;
	}
	locMap[iloc]= *it;
      }
      if(!bornPartner) {
	matched=false;
	break;
      }
      bornParticles.erase(bornPartner);
      // skip the next block if not enforcing colour consistency
      if(!enforceColourConsistency_) continue;
      if((**it).branchingParticle()->colourLine()) {
	if(cmap.find((**it).branchingParticle()->colourLine())!=cmap.end()) {
	  if(cmap[(**it).branchingParticle()->colourLine()]!=bornPartner->colourLine()) {
	    matched=false;
	  }
	}
	else {
	  cmap[(**it).branchingParticle()->colourLine()] = bornPartner->colourLine();
	}
      }
      if((**it).branchingParticle()->antiColourLine()) {
	if(cmap.find((**it).branchingParticle()->antiColourLine())!=cmap.end()) {
	  if(cmap[(**it).branchingParticle()->antiColourLine()]!=bornPartner->antiColourLine()) {
	    matched=false;
	  }
	}
	else {
	  cmap[(**it).branchingParticle()->antiColourLine()] = bornPartner->antiColourLine();
	}
      }
      // require a match
      if(!matched) break;
    }
    // if no match continue
    if(!matched) continue;
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
    if(forcePartners_) {
      locMap[emitter_  ]->colourPartner(locMap[spectator_]);
      locMap[spectator_]->colourPartner(locMap[emitter_  ]);
      locMap[emitter_  ]->branchingParticle()->partner(locMap[spectator_]->branchingParticle());
      locMap[spectator_]->branchingParticle()->partner(locMap[emitter_  ]->branchingParticle());
    }
    newTree.tree()->partnersSet(true);
    // set the beam particles
    PPair beams = lastXCombPtr()->lastParticles();
    // remove children of beams
    PVector beam_children = beams.first->children();
    if( (**newTree.tree()->incoming().begin()).branchingParticle()->momentum().z() /
    	beams.first->momentum().z() < 0.)
      swap( beams.first, beams.second );
    set<HardBranchingPtr>::iterator it = newTree.tree()->incoming().begin();
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
    // check the emitter and the spectator some how
    if(iemitter!=emitter_) continue;
    //do inverse momentum reconstruction
    if( !evolver()->showerModel()->kinematicsReconstructor()
	->deconstructHardJets( newTree.tree(), evolver(), ShowerInteraction::QCD ) ) continue;
    newTree.tree()->findNodes();
    newTree.weight(1.);
    hardTrees_.push_back( make_pair( newTree, 1. ) );
  }

  // select the tree
  PotentialTree chosen_hardTree;
  if (hardTrees_.size()==1) {
    chosen_hardTree = hardTrees_[0].first;
  }
  else {
    // if multiple trees pick the one with matching 
    // intermediate particle momenta
    for (unsigned int il=0; il<hardTrees_.size(); ++il){
      vector<pair <long int, Lorentz5Momentum> > particles;      
      PotentialTree testTree = hardTrees_[il].first;
      CKKWTreePtr check = testTree.tree();
      // get id and momenta of particles in hard tree
      for (set< HardBranchingPtr >::iterator it=check->branchings().begin();
	   it!=check->branchings().end(); ++it){
	particles.push_back(make_pair((*it)->branchingParticle()->id(),
				      (*it)->branchingParticle()->momentum()));
	if (!(*it)->children().empty()){
	  for (unsigned int ic=0; ic<(*it)->children().size(); ++ic)
	    particles.push_back(make_pair((*it)->branchingParticle()->children()[ic]->id(),
					  (*it)->branchingParticle()->
					  children()[ic]->momentum()));
	}	  	
	if ((*it)->parent()){
	  particles.push_back(make_pair((*it)->parent()->branchingParticle()->id(),
					(*it)->parent()->branchingParticle()->momentum()));
	  if (!(*it)->parent()->branchingParticle()->children().empty()){
	    for (unsigned int ic=0; 
		 ic<(*it)->parent()->branchingParticle()->children().size(); ++ic)
	      particles.push_back(make_pair((*it)->parent()->branchingParticle()->
					    children()[ic]->id(),
					    (*it)->parent()->branchingParticle()->
					    children()[ic]->momentum()));
	  }	    
	} 
      }

      // loop through and match to particles in real subprocess
      vector<pair <long int, Lorentz5Momentum> >::iterator part = particles.begin();
      // incoming
      for (; part!=particles.end(); ++part){
	if ((*part).first==real->incoming().first->id() &&
	    fuzzyEqual((*part).second, real->incoming().first->momentum()))
	  break;
      }
      if (part!=particles.end()) particles.erase(part);
      part = particles.begin();
      for (; part!=particles.end(); ++part){
	if ((*part).first==real->incoming().second->id() &&
	    fuzzyEqual((*part).second, real->incoming().second->momentum()))
	  break;
      }
      if (part!=particles.end()) particles.erase(part);
      // outgoing
      for (unsigned int io=0; io<real->outgoing().size(); ++io){
	part = particles.begin();
	for (; part!=particles.end(); ++part){
	  if ((*part).first==real->outgoing()[io]->id() &&
	      fuzzyEqual((*part).second, real->outgoing()[io]->momentum()))
	    break;
	}
	if (part!=particles.end()) particles.erase(part);
      }     
      // intermediate
      for (unsigned int ii=0; ii<real->intermediates().size(); ++ii){
	part = particles.begin();
	for (; part!=particles.end(); ++part){
	  if ((*part).first==real->intermediates()[ii]->id() &&
	      fuzzyEqual((*part).second, real->intermediates()[ii]->momentum()))
	    break;
	}
	if (part!=particles.end()) particles.erase(part);
      }     
      // intermediate CC with -1*momentum
      for (unsigned int ii=0; ii<real->intermediates().size(); ++ii){
	part = particles.begin();
	for (; part!=particles.end(); ++part){
	  if (!real->intermediates()[ii]->coloured() ||
	      (real->intermediates()[ii]->hasColour() &&
	       real->intermediates()[ii]->hasAntiColour())){
	    if ((*part).first==real->intermediates()[ii]->id() &&
		fuzzyEqual((*part).second, -1.*real->intermediates()[ii]->momentum()) )
	      break;
	  }
	  else {
	    if ((*part).first==-1.*real->intermediates()[ii]->id() &&
		fuzzyEqual((*part).second, -1.*real->intermediates()[ii]->momentum()) )
	      break;
	  }
	}
	if (part!=particles.end()) particles.erase(part);
      }
      // if all particles match, set as hardtree
      if (particles.empty()){
	chosen_hardTree = testTree;
	break;
      }     
    }
  }
  protoBranchings().clear();
  protoTrees().clear();
  hardTrees_.clear();
  if(! chosen_hardTree.tree() ) {
    return PotentialTree();
  }
  else
    return chosen_hardTree;
}

bool PowhegShowerHandler::checkDiagram(PotentialTree & tree,
				       tcDiagPtr loDiagram) const {
  set<HardBranchingPtr>::const_iterator cit;
  tcPDPair incoming;
  multiset<tcPDPtr,ParticleOrdering> outgoing;  
  //get the incoming and outgoing partons involved in hard process
  for( cit = tree.tree()->branchings().begin(); 
       cit != tree.tree()->branchings().end(); ++cit ){ 
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
  if(!incoming.first || !incoming.second)
    return 0.;
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
  if( tag.first  == loDiagram->getTag() || 
      tag.second == loDiagram->getTag() )
    tree.diagram(loDiagram);
  // check this is allowed
  return tree.diagram();
}


void PowhegShowerHandler::fillProtoTrees( ProtoTreePtr currentProtoTree,long id ) const {
  if(currentProtoTree->branchings().size()==(lastXCombPtr()->subProcess()->outgoing().size()+2)) 
    return;
  for( set<tProtoBranchingPtr>::const_iterator 
	 ita = currentProtoTree->branchings().begin();
       ita!=currentProtoTree->branchings().end();++ita) {
    for( set<tProtoBranchingPtr>::const_iterator 
	   itb = currentProtoTree->branchings().begin();
	 itb!=ita;++itb) {
      // can't merge two incoming branchings
      if( (**ita).status() == HardBranching::Incoming &&
	  (**itb).status() == HardBranching::Incoming ) continue;
      // if branching must be outgoing, skip incoming
      if(emitter_>=2 && ( (**ita).status() == HardBranching::Incoming || 
			 (**itb).status() == HardBranching::Incoming ))
	continue;
      // if branching must be incoming, skip outgoing
      if(emitter_<2 && ( (**ita).status() != HardBranching::Incoming &&
			(**itb).status() != HardBranching::Incoming ))
	continue;
      // get a new branching for this pair
      ProtoBranchingPtr currentBranching = getCluster(*ita,*itb);
      // check branching with the right PID
      if( ! currentBranching || 
      	  currentBranching->id() != id)
      	continue;
      // branching allowed so make a new Tree out of these branchings
      set< tProtoBranchingPtr > newTreeBranchings = currentProtoTree->branchings();
      newTreeBranchings.erase(*ita);
      newTreeBranchings.erase(*itb);
      newTreeBranchings.insert(currentBranching);
      ProtoTreePtr newProtoTree = new_ptr( ProtoTree( newTreeBranchings ) );
      // remove duplicate trees
      if( ! repeatProtoTree( newProtoTree ) ) protoTrees().insert( newProtoTree );
      // remove the current tree if it hasn't already been removed
      if( protoTrees().find( currentProtoTree ) != protoTrees().end() )
	protoTrees().erase( currentProtoTree );
      // do recursion
      fillProtoTrees( newProtoTree , id);
    }
  }
}

bool  PowhegShowerHandler::repeatProtoTree( ProtoTreePtr currentProtoTree ) const {
  // loop over all prototrees and see 
  // how many ProtoBranchings of current ProtoTree are found in each
  for( set< ProtoTreePtr >::const_iterator cit = protoTrees().begin();
       cit != protoTrees().end(); ++cit ) {
    unsigned int no_matches = 0;
    for( set< tProtoBranchingPtr >::const_iterator ckt 
	   = currentProtoTree->branchings().begin(); 
	 ckt != currentProtoTree->branchings().end(); ckt++ ) {
      if( (*cit)->branchings().find( *ckt ) != (*cit)->branchings().end() )
	++no_matches;
    }
    // return true if all match
    if( no_matches == currentProtoTree->branchings().size() )
      return true;
  }
  return false;
}

tProtoBranchingPtr PowhegShowerHandler::getCluster( tProtoBranchingPtr b1,
						    tProtoBranchingPtr b2 ) const {
  // look for the clustered pair in protoBranchings_
  for(set<ProtoBranchingPtr>::const_iterator cit = protoBranchings().begin();
      cit != protoBranchings().end(); ++cit) {
    // both outgoing
    if(b1->status()==HardBranching::Outgoing &&
       b2->status()==HardBranching::Outgoing) {
      if((**cit).status()!=HardBranching::Outgoing||
	 (**cit).children().empty()) continue;
      if( ( b1 == (**cit).children()[0] && b2 == (**cit).children()[1] ) ||
	  ( b1 == (**cit).children()[1] && b2 == (**cit).children()[0] ) )
	return *cit;
    }
    // first incoming
    else if(b1->status()==HardBranching::Incoming) {
      if((**cit).backChildren().empty() ) continue;
      if(b1!=(**cit).backChildren()[0]) continue;
      if(b2==(**cit).backChildren()[1]) return *cit;
    }
    // second incoming
    else if(b2->status()==HardBranching::Incoming) {
      if((**cit).backChildren().empty() ) continue;
      if(b2!=(**cit).backChildren()[0]) continue;
      if(b1==(**cit).backChildren()[1]) return *cit;
    }
  }
  // is branching incoming or outgoing
  bool incoming = b1->status()==HardBranching::Incoming ||
                  b2->status()==HardBranching::Incoming;
  // get the branching
  BranchingElement theBranching;
  if( !incoming ) theBranching =   allowedFinalStateBranching( b1, b2 );
  else            theBranching = allowedInitialStateBranching( b1, b2 );

  //if branching is not allowed return null ProtoBrancing
  if( !theBranching.first )
    return ProtoBranchingPtr();

  // get the ParticleData object for the new branching
  tcPDPtr particle_data = incoming ?
    getParticleData( theBranching.second[1] ) : getParticleData( theBranching.second[0] );

  // create clustered ProtoBranching
  ProtoBranchingPtr clusteredBranch;

  // outgoing
  if( !incoming ){
    Lorentz5Momentum pairMomentum = b1->momentum() + b2->momentum(); 
    pairMomentum.setMass(ZERO);
    clusteredBranch = new_ptr(ProtoBranching(particle_data,HardBranching::Outgoing,
					     pairMomentum, theBranching.first));
    if(particle_data->iColour()==PDT::Colour0)
      return ProtoBranchingPtr();
    else if(particle_data->iColour()==PDT::Colour3) {
      if(b1->particle()->iColour()==PDT::Colour3 && b2->particle()->iColour()==PDT::Colour8) {
	if(b1->colourLine()!=b2->antiColourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->colourLine(b2->colourLine());
      }
      else if(b2->particle()->iColour()==PDT::Colour3 && b1->particle()->iColour()==PDT::Colour8) {
	if(b2->colourLine()!=b1->antiColourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->antiColourLine(b1->colourLine());
      }
      else
	assert(false);
      clusteredBranch->type(ShowerPartnerType::QCDColourLine);
    }
    else if(particle_data->iColour()==PDT::Colour3bar) {
      if(b1->particle()->iColour()==PDT::Colour3bar && b2->particle()->iColour()==PDT::Colour8) {
	if(b1->antiColourLine()!=b2->colourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->antiColourLine(b2->antiColourLine());
      }
      else if(b2->particle()->iColour()==PDT::Colour3bar && b1->particle()->iColour()==PDT::Colour8) {
	if(b2->antiColourLine()!=b1->colourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->antiColourLine(b1->antiColourLine());
      }
      else
	assert(false);
      clusteredBranch->type(ShowerPartnerType::QCDAntiColourLine);
    }
    else if(particle_data->iColour()==PDT::Colour8) {
      tProtoBranchingPtr coloured,antiColoured;
      if(b1->particle()->iColour()==PDT::Colour3 &&
	 b2->particle()->iColour()==PDT::Colour3bar) {
	coloured     = b1;
	antiColoured = b2;
      }
      else if(b2->particle()->iColour()==PDT::Colour3 &&
	      b1->particle()->iColour()==PDT::Colour3bar) {
	coloured     = b2;
	antiColoured = b1;
      }
      else if(b1->particle()->iColour()==PDT::Colour8 &&
	      b2->particle()->iColour()==PDT::Colour8 ) {
	if(b1->colourLine()==b2->antiColourLine()) {
	  coloured     = b2;
	  antiColoured = b1;
	}
	else if(b2->colourLine()==b1->antiColourLine()) {
	  coloured     = b1;
	  antiColoured = b2;
	}
	else
	  return ProtoBranchingPtr();
      } 
      else
	assert(false);
      // softest particle is the emitted
      if(coloured->momentum().t()>antiColoured->momentum().t()) {
	clusteredBranch->type(ShowerPartnerType::QCDAntiColourLine);
      }
      else {
	clusteredBranch->type(ShowerPartnerType::QCDColourLine);
      }
    }
    else
      assert(false);
  }
  // incoming
  else {
    Lorentz5Momentum pairMomentum = b1->momentum() - b2->momentum();
    pairMomentum.setMass( ZERO );
    // check for CC
    if( particle_data->CC() &&
	( b1->id() != theBranching.second[0] ||
	  b2->id() != theBranching.second[2] ) ) {
      particle_data = particle_data->CC();
    }
    clusteredBranch = new_ptr(ProtoBranching(particle_data,HardBranching::Incoming,
					     pairMomentum,theBranching.first));
    // work out the type of branching
    if(b1->particle()->iColour()==PDT::Colour3) {
      b1->type(ShowerPartnerType::QCDColourLine);
      if(b2->particle()->iColour()==PDT::Colour3 &&
	 particle_data->iColour()==PDT::Colour8) {
	clusteredBranch->    colourLine(b1->colourLine());
	clusteredBranch->antiColourLine(b2->colourLine());
      }
      else if(b2->particle()->iColour()==PDT::Colour8 &&
	      particle_data->iColour()==PDT::Colour3) {
	if(b1->colourLine()!=b2->colourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->colourLine(b2->antiColourLine());
      }
      else
	assert(false);
    }
    else if(b1->particle()->iColour()==PDT::Colour3bar) {
      b1->type(ShowerPartnerType::QCDAntiColourLine);
      if(b2->particle()->iColour()==PDT::Colour3bar &&
	 particle_data->iColour()==PDT::Colour8) {
	clusteredBranch->    colourLine(b2->antiColourLine());
	clusteredBranch->antiColourLine(b1->antiColourLine());
      }
      else if(b2->particle()->iColour()==PDT::Colour8 &&
	      particle_data->iColour()==PDT::Colour3bar) {
	if(b1->antiColourLine()!=b2->antiColourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->antiColourLine(b2->colourLine());
      }
      else
	assert(false);
    }
    else if(b1->particle()->iColour()==PDT::Colour8) {
      if(b2->particle()->iColour()==PDT::Colour3) {
	if(b1->colourLine()!=b2->colourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->antiColourLine(b1->antiColourLine());	
	b1->type(ShowerPartnerType::QCDColourLine);
      }
      else if(b2->particle()->iColour()==PDT::Colour3bar) {
	if(b1->antiColourLine()!=b2->antiColourLine())
	  return ProtoBranchingPtr();
	clusteredBranch->    colourLine(b1->colourLine());
	b1->type(ShowerPartnerType::QCDAntiColourLine);
      }
      else if(b2->particle()->iColour()==PDT::Colour8) {
	if(b1->colourLine()==b2->colourLine()) {	
	  b1->type(ShowerPartnerType::QCDColourLine);
	  clusteredBranch->antiColourLine(b1->antiColourLine());
	  clusteredBranch->colourLine(b2->antiColourLine());
	}
	else if(b1->antiColourLine()==b2->antiColourLine()) {
	  b1->type(ShowerPartnerType::QCDAntiColourLine);
	  clusteredBranch->    colourLine(b1->colourLine());
	  clusteredBranch->antiColourLine(b2->colourLine());
	}
	else {
	  return ProtoBranchingPtr();
	}
      }
      else
	assert(false);
    }
    else
      assert(false);
  }
  protoBranchings().insert(clusteredBranch);
  //set children relations 
  // outgoing
  if( !incoming ){
    clusteredBranch->addChild( b1 );	    
    clusteredBranch->addChild( b2 );  
  }
  else {
    clusteredBranch->addBackChild( b1 );	    
    clusteredBranch->addBackChild( b2 );
  }
  return clusteredBranch;
}


BranchingElement PowhegShowerHandler::
allowedFinalStateBranching( tProtoBranchingPtr & b1, tProtoBranchingPtr & b2) const {
  // check with normal ID's
  pair< long, long > ptest = make_pair( b1->id(), b2->id() );
  map< pair< long, long >, pair< SudakovPtr, IdList > >::const_iterator 
    split = allowedFinal_.find(ptest);
  if( split != allowedFinal_.end() ) {
    if(  split->second.second[1] != ptest.first ) swap( b1, b2 );
    return split->second;
  }
  // check with CC
  if( b1->particle()->CC() ) ptest.first  *= -1;
  if( b2->particle()->CC() ) ptest.second *= -1;
  split = allowedFinal_.find( ptest );
  if( split != allowedFinal_.end() ) {
    // cc the idlist only be for qbar g clusterings
    BranchingElement ccBranch = split->second;
    if( getParticleData( ccBranch.second[0] )->CC() ) ccBranch.second[0] *= -1;
    if( getParticleData( ccBranch.second[1] )->CC() ) ccBranch.second[1] *= -1;
    if( getParticleData( ccBranch.second[2] )->CC() ) ccBranch.second[2] *= -1;
    if( split->second.second[1] !=  ptest.first ) swap( b1, b2);
    return ccBranch;
  }
  // not found found null pointer
  return make_pair( SudakovPtr(), IdList() );
}

BranchingElement PowhegShowerHandler::allowedInitialStateBranching( tProtoBranchingPtr & b1,
								    tProtoBranchingPtr & b2) const {
  if(b2->status()==HardBranching::Incoming) swap(b1,b2);
  // is initial parton an antiparticle
  bool cc = b1->id() < 0;
  //gives range of allowedInitial_ with matching first abs( id )
  pair< multimap< long, pair< SudakovPtr, IdList > >::const_iterator,
    multimap< long, pair< SudakovPtr, IdList > >::const_iterator >
    location = allowedInitial_.equal_range( abs( b1->id() ) );
  //iterates over this range
  for( multimap< long, pair< SudakovPtr, IdList> >::const_iterator it = location.first;
       it != location.second; ++it ) {
    //test id for second particle in pair
    long idtest = it->second.second[2];
    //if it is antiparticle *= -1
    if( cc && getParticleData( idtest )->CC() ) idtest *= -1;
    // does second id match the test
    if( idtest == b2->id() ) return it->second;
    //if the the IS parton is a gluon and charge conjugate of second parton mathes accept
    if( idtest == -b2->id() &&
        ! b1->particle()->CC() ) return it->second;
  }
  // not found found null pointer
  return make_pair(SudakovPtr(),IdList());
}

bool PowhegShowerHandler::fuzzyEqual(Lorentz5Momentum  a, 
				     Lorentz5Momentum  b) const{
  // check momenta are within 1% of each other
  if ( (a.e()==ZERO && b.e()==ZERO) || (a.e()/b.e()>0.99 && a.e()/b.e()<1.01) ){
    if ((a.x()==ZERO && b.x()==ZERO) || (a.x()/b.x()>0.99 && a.x()/b.x()<1.01) ){
      if ((a.y()==ZERO && b.y()==ZERO) || (a.y()/b.y()>0.99 && a.y()/b.y()<1.01) ){
	if ((a.z()==ZERO && b.z()==ZERO) || (a.z()/b.z()>0.99 && a.z()/b.z()<1.01) )
	  return true;
      }
    }
  }
  return false;
}

void PowhegShowerHandler::doinit() {
  ShowerHandler::doinit();
  // extract the allowed branchings
  // final-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->finalStateBranchings().begin();
      it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    pair<long,long> prod(make_pair(it->second.second[1],it->second.second[2]));
    allowedFinal_.insert(make_pair(prod,it->second));
    swap(prod.first,prod.second);
    allowedFinal_.insert(make_pair(prod,it->second));
  }
  // initial-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->initialStateBranchings().begin();
      it != evolver()->splittingGenerator()->initialStateBranchings().end(); ++it) {
    allowedInitial_.insert(make_pair(it->second.second[0],it->second));
  }

}

void PowhegShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << theFactory << allowedInitial_ << allowedFinal_
     << subtractionIntegral_ << enforceColourConsistency_ << forcePartners_;
}

void PowhegShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> theFactory >> allowedInitial_ >> allowedFinal_
     >> subtractionIntegral_ >> enforceColourConsistency_ >> forcePartners_;
}


// Static variable needed for the type description system in ThePEG.
DescribeClass<PowhegShowerHandler,Herwig::ShowerHandler>
describeHerwigPowhegShowerHandler("Herwig::PowhegShowerHandler", 
				  "HwMatchbox.so HwMatching.so");

void PowhegShowerHandler::Init() {

  static ClassDocumentation<PowhegShowerHandler> documentation
    ("The PowhegShowerHandler class");

  static Reference<PowhegShowerHandler,MatchboxFactory> interfaceFactory
    ("Factory",
     "The factory object to use.",
     &PowhegShowerHandler::theFactory, false, false, true, false, false);

  static Switch<PowhegShowerHandler,bool> interfaceEnforceColourConsistency
    ("EnforceColourConsistency",
     "Force the Born and real emission colour flows to be consistent",
     &PowhegShowerHandler::enforceColourConsistency_, false, false, false);
  static SwitchOption interfaceEnforceColourConsistencyYes
    (interfaceEnforceColourConsistency,
     "Yes",
     "Enforce the consistency",
     true);
  static SwitchOption interfaceEnforceColourConsistencyNo
    (interfaceEnforceColourConsistency,
     "No",
     "Don't enforce consistency",
     false);

  static Switch<PowhegShowerHandler,bool> interfaceForcePartners
    ("ForcePartners",
     "Whether or not to force the partners to be those from the kinematic generation",
     &PowhegShowerHandler::forcePartners_, false, false, false);
  static SwitchOption interfaceForcePartnersYes
    (interfaceForcePartners,
     "Yes",
     "Force them",
     true);
  static SwitchOption interfaceForcePartnersNo
    (interfaceForcePartners,
     "No",
     "Don't force them",
     false);

}
