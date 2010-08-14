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
   * Fix the fwd branching connections for spacelike shower
   */
  bool fixFwdBranchings();

  /**
   * Access the incoming branchings
   */
  set<HardBranchingPtr> & incoming() {return _spacelike;}

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type interaction() {return _interaction;}

  /**
   *  Get LowestPt in which ever jet definition from the shower variables
   */
  Energy lowestPt( int jetMeasureMode, Energy2 s );

  /**
   *  Get lowest Pt in which ever jet definition from the hardtree momentum
   */
  Energy lowestPtMomentum( int jetMeasureMode, int cutOption );
  
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
   * Access the intermediate lines
   */
  map< HardBranchingPtr, HardBranchingPtr > & getIntermediates() { return _theIntermediates; }



  /**
   * Returns true if all lines in tree are ordered in /tilde(q)
   */
  bool checkHardOrdering();
  
  /**
   * Returns true if all spacelike lines are ordered in x
   */
  bool checkXOrdering();

  /**
   * Calls recursive function to fill externals and nodes
   * then finds the internal lines from the nodes
   */
  bool findNodes( );

  /**
   * Returns sum of pts of all branchings in tree
   */
  Energy totalPt() { return _total_pt; } 

  LorentzRotation showerRot() { return _showerRot; }

  void showerRot( LorentzRotation r ) { _showerRot = r; }

  /**
   * Remove all back child relations
   * Needs to be done once a hardTree has been selected or rejected 
   * (i.e. once back child relations are no longer required) so that
   * memory leaks associated with cyclic pointer relations are avoided.
   */
  void removeBackChildren();

  /**
   *  Whether or not the evolution partners are set
   */
  bool partnersSet() const {return partnersSet_;}

  /**
   *  Whether or not the evolution partners are set
   */
  void partnersSet(bool in) {partnersSet_=in;}

private:

  /**
   * Recursive function to fix the parent assignments
   */
  bool fixParents( HardBranchingPtr );

  /**
   * Recursive function to fill externals, nodes and intermediates from the time-like showers
   */
  bool fillNodes( HardBranchingPtr );

  /**
   * Function to recursively find the hard line scales
   **/
  void fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine );

private:

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type _interaction;


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
   * Map of the start and the end branchings of an intermediate line
   */
  map< HardBranchingPtr, HardBranchingPtr > _theIntermediates;

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
   * The lowest pt of the branchings in the hardtree in whatever
   * jet measure according to the hardtree momenta (not the shower variables)
   */
  Energy  _lowestPtMomentum;

  /**
   *  The sum of the pts of all branchings
   */
  Energy _total_pt;

  //rotation to shower frame
  LorentzRotation _showerRot;

  /**
   *  Whether or not partners are set
   */
  bool partnersSet_;

};

  /**
   *  Output operator for testing
   */
  ostream & operator << (ostream &, const HardTree & );

}

#endif /* HERWIG_HardTree_H */
