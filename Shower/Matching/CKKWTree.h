// -*- C++ -*-
#ifndef HERWIG_CKKWTree_H
#define HERWIG_CKKWTree_H
//
// This is the declaration of the CKKWTree class.
//

#include "Herwig/Shower/Base/HardTree.h"
#include "CKKWTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the CKKWTree class.
 */
class CKKWTree: public HardTree {

public:

  /**
   * The default constructor.
   */
  CKKWTree(vector<HardBranchingPtr>, vector<HardBranchingPtr>, ShowerInteraction::Type);

public:

  /**
   * Returns true if all spacelike lines are ordered in x
   */
  bool checkXOrdering();

  /**
   *  Get lowest Pt in which ever jet definition from the hardtree momentum
   */
  Energy lowestPtMomentum( int jetMeasureMode, int cutOption );

  /**
   * Returns true if all lines in tree are ordered in /tilde(q)
   */
  bool checkHardOrdering();

  /**
   * Returns sum of pts of all branchings in tree
   */
  Energy totalPt() { return totalpT_; } 

  /**
   * Calls recursive function to fill externals and nodes
   * then finds the internal lines from the nodes
   */
  void findNodes();

  /**
   * Access the nodal branchings
   */
  map< HardBranchingPtr, Energy > & getNodes() { return nodes_; }

  /**
   *  Get LowestPt in which ever jet definition from the shower variables
   */
  Energy lowestPt( int jetMeasureMode, Energy2 s );

  /**
   *  Is the tree ordered?
   */
  bool ordered() const {return ordered_;}

protected:

  /**
   * Recursive function to find the lowest jet measure in a hardtree from clustered momenta
   **/
  void getLowestJetMeasure( HardBranchingPtr branch, int jetMeasureMode, int cutOption );

  /**
   * Function for finding the hadronic jet measure of two partons
   **/
  Energy hadronJetMeasure( const Lorentz5Momentum & p1, 
			   const Lorentz5Momentum & p2,
			   bool final );
  /**
   * Function for finding the Durham or LUCLUS jet measures of two partons
   **/
  Energy getJetMeasure( const Lorentz5Momentum & p1,
			const Lorentz5Momentum & p2,
			int jetMeasureMode );

  /**
   * Function to determine whether a branching consists of external partons
   **/
  bool externalBranching( HardBranchingPtr a, HardBranchingPtr b );

  /**
   * Function to recursively find the hard line scales
   **/
  void fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine );

  /**
   * Recursive function to fill externals, nodes and intermediates from the time-like showers
   */
  bool fillNodes( HardBranchingPtr );

  /**
   * Recursive function to fix the parent assignments
   */
  bool fixParents( HardBranchingPtr );

private:

  /**
   * Map containing all nodes with the ingoing partons and their scale 
   * (this is the the intermediates and their ending node).
   */
  map< HardBranchingPtr,  Energy > nodes_;

  /**
   * Scales and z along each hard line to check ordering
   */
  vector< vector< pair< Energy, double > > > hardLineScales_;

  /**
   * The lowest pt of the branchings in the hardtree in whatever
   * jet measure according to the hardtree momenta (not the shower variables)
   */
  Energy  lowestpTMomentum_;

  /**
   *  The sum of the pts of all branchings
   */
  Energy totalpT_;

  /**
   *  The hardBranching of softest branching
   *  This is found by looking at tree end points in fillNodes
   */
  HardBranchingPtr  lowestpT_;

  /**
   *  Is the tree ordered
   */
  bool ordered_;
  
};

}

#endif /* HERWIG_CKKWTree_H */
