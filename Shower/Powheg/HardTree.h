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
  friend ostream & operator<<(ostream &, const HardTree & );

public:

  /**
   * The default constructor.
   */
  HardTree(vector<HardBranchingPtr>,vector<HardBranchingPtr>);

  /**
   *  Match particles in the ShowerTree to branchings in the HardTree
   */
  inline void connect(ShowerParticlePtr,HardBranchingPtr);
  
  /**
   *  Access the map between the ShowerParticle and the HardBranching
   */
  inline map<ShowerParticlePtr,tHardBranchingPtr> & particles();

  /**
   *  Access the set of branchings
   */
  inline set<HardBranchingPtr> & branchings();
  
  /**
   * Access the incoming branchings
   */
  inline set<HardBranchingPtr> & incoming();

private:

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
};

/**
 * The HardBranching class is designed to contain the information needed for
 * an individual branching in the POWHEG approach
 */
class HardBranching : public Base {

  /**
   *  The HardTree is friend
   */
  friend class HardTree;

public:

  /**
   * The default constructor
   * @param particle The particle which is branching
   * @param sudakov  The Sudakov form factor for the branching
   * @param parent   The parent for the branching
   * @param incoming Whether the particle is incoming or outgoing
   */
  inline HardBranching(ShowerParticlePtr particle,
			SudakovPtr sudakov,
			tHardBranchingPtr parent,bool incoming);

  /**
   * Add a child of the branching
   * @param child The child of the branching
   */
  inline void addChild(HardBranchingPtr child);

  /**
   * Return the ShowerParticlePtr of the branching particle.
   */
  inline ShowerParticlePtr branchingParticle();

  void setMomenta(LorentzRotation R,double alpha,Lorentz5Momentum pt);

  /**
   *  The branching particle
   */
  ShowerParticlePtr _particle;

  /**
   *  The rescaled momentum
   */
  Lorentz5Momentum _original;

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum _p;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum _n;

  /**
   *  The transverse momentum vector
   */
  Lorentz5Momentum _qt;

  /**
   *  The momentum the particle should have as the start of a shower
   */
  Lorentz5Momentum _shower;

  /**
   *  The transverse momentum
   */
  Energy _pt;

  /**
   *  Whether the branching is incoming or outgoing
   */
  bool _incoming;

  /**
   *  Information on the Shower variables for the branching
   */
  //@{
  /**
   *  The evolution scale
   */
  Energy _scale;

  /**
   *  The energy fraction
   */
  double _z;

  /**
   *  The azimthual angle
   */
  double _phi;
  //@}

  /**
   *  The parent of the branching
   */
  tHardBranchingPtr _parent;

  /**
   *  The Sudakov form-factor
   */
  SudakovPtr _sudakov;

  /**
   * The children
   */
  vector<HardBranchingPtr> _children;

  /**
   *  The beam particle
   */
  PPtr _beam;
};

inline ostream & operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x._branchings.begin();
      it!=x._branchings.end();++it) {
    os << "Hard Particle: " << *(**it)._particle << " has Sudakov " 
       << (**it)._sudakov << "\n";
    os << "It's colour lines are " << (**it)._particle->colourLine() << "\t" 
       <<  (**it)._particle->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it)._children.size();++iy) {
      os << "\t Children : " << *(**it)._children[iy]->_particle
	 << "\n";
      os << "It's colour lines are " << (**it)._children[iy]->_particle->colourLine() << "\t" 
	 <<  (**it)._children[iy]->_particle->antiColourLine() << "\n";
    }
  }
  for(set<HardBranchingPtr>::const_iterator it=x._spacelike.begin();
      it!=x._spacelike.end();++it) {
    os << "SpaceLike: " << *(**it)._particle << " has Sudakov" 
       << (**it)._sudakov << "\n";
    for(unsigned int iy=0;iy<(**it)._children.size();++iy) {
      os << "\t Children: " << *(**it)._children[iy]->_particle
	 << "\n";
      os << "It's colour lines are " << (**it)._children[iy]->_particle->colourLine() << "\t" 
	 <<  (**it)._children[iy]->_particle->antiColourLine() << "\n";
    }
  }
  return os;
}

}

#include "HardTree.icc"

#endif /* HERWIG_HardTree_H */
