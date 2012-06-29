// -*- C++ -*-
#ifndef HERWIG_HardTree_H
#define HERWIG_HardTree_H
//
// This is the declaration of the HardTree class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "HardBranching.h"
#include "HardTree.fh"

namespace Herwig {

using namespace ThePEG;
/**
 * The HardTree class is designed to contain the information required
 * to implement the POWHEG approach for Monte Carlo at next-to-leading order.
 */
class HardTree : public Base {

  /**
   *  Output operator for testing
   */
  friend ostream & operator << (ostream &, const HardTree & );

public:

  /**
   * The default constructor.
   */
  HardTree(vector<HardBranchingPtr>, vector<HardBranchingPtr>, ShowerInteraction::Type);

  /**
   *  Match particles in the ShowerTree to branchings in the HardTree
   */
  bool connect(ShowerParticlePtr particle, HardBranchingPtr branching) {
    if(branchings_.find(branching)==branchings_.end()) return false;
    particles_[particle]=branching;
    return true;
  }

  /**
   *  Match the prticles in the ShowerTree to the branchings in the HardTree
   */
  bool connect(ShowerTreePtr);
  
  /**
   *  Access the map between the ShowerParticle and the HardBranching
   */
  map<ShowerParticlePtr,tHardBranchingPtr> & particles() 
  {return particles_;}

  /**
   *  Access the set of branchings
   */
  set<HardBranchingPtr> & branchings() {return branchings_;}

  /**
   * Access the incoming branchings
   */
  set<HardBranchingPtr> & incoming() {return spacelike_;}

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type interaction() {return interaction_;}

  /**
   *  Get the Rotation to be applied to the tree
   */
  LorentzRotation showerRot() { return showerRot_; }

  /**
   *  Set the Rotation to be applied to the tree
   */
  void showerRot( LorentzRotation r ) { showerRot_ = r; }

private:

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type interaction_;

  /**
   *  The ShowerTree
   */
  ShowerTreePtr _tree;

  /**
   *  Map from the particles in the ShowerTree to the HardBranchings
   */
  map<ShowerParticlePtr,tHardBranchingPtr> particles_;

  /**
   *  The HardBranchings in the hard process
   */
  set<HardBranchingPtr> branchings_;

  /**
   *  The HardBranchings which initiate the space-like showers
   */
  set<HardBranchingPtr> spacelike_;

  /**
   * Rotation to shower frame
   */
  LorentzRotation showerRot_;

};

  /**
   *  Output operator for testing
   */
  ostream & operator << (ostream &, const HardTree & );

}

#endif /* HERWIG_HardTree_H */
