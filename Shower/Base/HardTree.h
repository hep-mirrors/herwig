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
  HardTree(vector<HardBranchingPtr>,vector<HardBranchingPtr>);

  /**
   *  Match particles in the ShowerTree to branchings in the HardTree
   */
  bool connect(ShowerParticlePtr particle, HardBranchingPtr branching) {
    if(_branchings.find(branching)==_branchings.end()) return false;
    _particles[particle]=branching;
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
  {return _particles;}

  /**
   *  Access the set of branchings
   */
  set<HardBranchingPtr> & branchings() {return _branchings;}
  
  /**
   * Access the incoming branchings
   */
  set<HardBranchingPtr> & incoming() {return _spacelike;}

  /**
   *  Get LowestPt in which ever jet definition
   */
  Energy lowestPt( int jetMeasureMode, Energy2 s );
  
  /**
   * Access the external branchings
   */
  map< ShowerParticlePtr, HardBranchingPtr > & getExternals() 
  { return _theExternals; }

  /**
   * Access the nodal branchings
   */
  map< HardBranchingPtr, Energy > & getNodes() { return _theNodes; }

  /**
   * Access the internal lines
   */
  map< long, pair< Energy, Energy > > & getInternals() { return _theInternals; }



  /**
   * Returns true if all lines in tree are ordered in /tilde(q)
   */
  bool checkHardOrdering();
  
  /**
   * Calls recursive function to fill externals and nodes
   * then finds the internal lines from the nodes
   */
  bool findNodes( );

  /**
   * Returns sum of pts of all branchings in tree
   */
  Energy totalPt() { return _total_pt; } 


private:

 /**
   * Recursive function to fill externals and nodes, also connects the parents
   */
  bool fillNodes( HardBranchingPtr, HardBranchingPtr );

  /**
   * Function to recursively find the hard line scales
   **/
  void fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine );

  /**
   * Scales and z along each hard line to check ordering
   */
  vector< vector< pair< Energy, double > > > _hard_line_scales;

  /**
   *  The ShowerTree
   */
  ShowerTreePtr _tree;

  /**
   *  Map from the particles in the ShowerTree to the HardBranchings
   */
  map<ShowerParticlePtr,tHardBranchingPtr> _particles;

  /**
   *  The HardBranchings in the hard process
   */
  set<HardBranchingPtr> _branchings;

  /**
   *  The HardBranchings which initiate the space-like showers
   */
  set<HardBranchingPtr> _spacelike;

  /**
   * Map containing external particles and their branchings
   */
  map< ShowerParticlePtr, HardBranchingPtr > _theExternals;

  /**
   * Map containing all nodes with the ingoing partons and their scale 
   * (this is the the intermediates and their ending node).
   */
  map< HardBranchingPtr,  Energy > _theNodes;
  
  /**
   * Map containing all internal line ids and 
   * their start and end node scale, qtilde.
   */
  map< long, pair< Energy, Energy > > _theInternals;

  /**
   *  The hardBranching of softest branching
   *  This is found by looking at tree end points in fillNodes
   */
   HardBranchingPtr  _lowestPt;

  /**
   *  The sum of the pts of all branchings
   */
  Energy _total_pt;

};


/**
 *  Output operator for testing
 */
ostream & operator << (ostream &, const HardTree & );

}

#endif /* HERWIG_HardTree_H */
