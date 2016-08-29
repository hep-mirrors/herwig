// -*- C++ -*-
#ifndef HERWIG_ProtoTree_H
#define HERWIG_ProtoTree_H
//
// This is the declaration of the ProtoTree class.
//

#include "ProtoBranching.h"
#include "CKKWTree.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Declare the pointers
 */
class ProtoTree;
ThePEG_DECLARE_POINTERS(Herwig::ProtoTree,ProtoTreePtr);

/**
 *  Class for a prototype tree
 */
class ProtoTree:public Base {

public:

  /**
   * Default constructor
   */  
  ProtoTree() {}
  
  /**
   *  Constructor
   */
  ProtoTree(const set<tProtoBranchingPtr> & newBranchings) :
    branchings_( newBranchings ) {}
  
  /**
   *  Add a branching
   */
  void addBranching( tProtoBranchingPtr Branching ){
    branchings_.insert( Branching );
  }
  
  const set< tProtoBranchingPtr > & branchings() const {
    return  branchings_;
  }

  /**
   *  Create the HardTree
   */
  CKKWTreePtr createHardTree() {
    vector<HardBranchingPtr> branchings,spacelike;
    map<ColinePtr,ColinePtr> cmap;
    for(set<tProtoBranchingPtr>::const_iterator it=branchings_.begin();
	it!=branchings_.end();++it) {
      if((**it).status()==HardBranching::Outgoing) {
	branchings.push_back(createTimeLikeBranching(*it,cmap));
      }
      else {
	HardBranchingPtr space;
	branchings.push_back(createSpaceLikeBranching(*it,space,cmap));
	spacelike.push_back(space);
      }
    }
    // Create the hard tree
    return new_ptr(CKKWTree( branchings, spacelike , ShowerInteraction::QCD));
  }

protected:

  /**
   *  Create a timelike branching
   */
  HardBranchingPtr createTimeLikeBranching(tProtoBranchingPtr branch,map<ColinePtr,ColinePtr> & cmap) {
    ShowerParticlePtr particle = new_ptr( ShowerParticle( branch->particle() , true ) );
    particle->set5Momentum( branch->momentum() );
    HardBranchingPtr newBranch = new_ptr( HardBranching( particle, branch->sudakov(),
							 HardBranchingPtr(), 
							 HardBranching::Outgoing ) );
    if(branch->colourLine()) {
      if(cmap.find((branch->colourLine()))==cmap.end())
	cmap[branch->colourLine()] = new_ptr(ColourLine());
      cmap[branch->colourLine()]->addColoured(particle);
    }
    if(branch->antiColourLine()) {
      if(cmap.find((branch->antiColourLine()))==cmap.end())
	cmap[branch->antiColourLine()] = new_ptr(ColourLine());
      cmap[branch->antiColourLine()]->addAntiColoured(particle);
    }
    Lorentz5Momentum pnew;
    if(branch->children().empty()) {
      pnew = branch->momentum();
    }
    else {
      for(unsigned int ix=0;ix<branch->children().size();++ix) {
	HardBranchingPtr child = createTimeLikeBranching(branch->children()[ix],cmap);
	newBranch->addChild(child);
	child->parent(newBranch);
	pnew += child->branchingParticle()->momentum();
	pnew.setMass(branch->particle()->mass());
      }
    }
    particle->set5Momentum( pnew );
    if(branch->type()!=ShowerPartnerType::Undefined)
      newBranch->type(branch->type());
    return newBranch;
  }

  /**
   *  Create a spacelike branching
   */
  HardBranchingPtr createSpaceLikeBranching(tProtoBranchingPtr branch,
					    HardBranchingPtr & spacelike,
					    map<ColinePtr,ColinePtr> & cmap) {
    ShowerParticlePtr particle = new_ptr( ShowerParticle( branch->particle() , false ) );
    particle->set5Momentum( branch->momentum() );
    if(branch->colourLine()) {
      if(cmap.find((branch->colourLine()))==cmap.end())
	cmap[branch->colourLine()] = new_ptr(ColourLine());
      cmap[branch->colourLine()]->addColoured(particle);
    }
    if(branch->antiColourLine()) {
      if(cmap.find((branch->antiColourLine()))==cmap.end())
	cmap[branch->antiColourLine()] = new_ptr(ColourLine());
      cmap[branch->antiColourLine()]->addAntiColoured(particle);
    }
    HardBranchingPtr newBranch = new_ptr( HardBranching( particle, SudakovPtr(),
							 HardBranchingPtr(), 
							 HardBranching::Incoming ) );
    if(!branch->backChildren().empty()) {
      HardBranchingPtr newTimeLike  =  createTimeLikeBranching(branch->backChildren()[1],cmap);
      HardBranchingPtr newSpaceLike = createSpaceLikeBranching(branch->backChildren()[0],
							       spacelike,cmap);
      newBranch  ->parent(newSpaceLike);
      newTimeLike->parent(newSpaceLike);
      newSpaceLike->addChild(newBranch);
      newSpaceLike->addChild(newTimeLike);
      newSpaceLike->sudakov(branch->sudakov());
    }
    else {
      spacelike = newBranch;
    }
    if(branch->type()!=ShowerPartnerType::Undefined)
      newBranch->type(branch->type());
    return newBranch;
  }


private:
  
  /**
   *  The branchings in the tree
   */
  set< tProtoBranchingPtr > branchings_;

};

}

#endif /* HERWIG_ProtoTree_H */
