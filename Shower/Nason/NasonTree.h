// -*- C++ -*-
#ifndef HERWIG_NasonTree_H
#define HERWIG_NasonTree_H
//
// This is the declaration of the NasonTree class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "NasonTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The NasonTree class is designed to contain the information required
 * to implement the Nason approach for Monte Carlo at next-to-leading order.
 */
class NasonTree : public Base {

public:

  /**
   * The default constructor.
   */
  NasonTree(vector<NasonBranchingPtr>,vector<NasonBranchingPtr>);

  /**
   *  Match particles in the ShowerTree to branchings in the NasonTree
   */
  inline void connect(ShowerParticlePtr,NasonBranchingPtr);
  
  /**
   *  Access the map between the ShowerParticle and the NasonBranching
   */
  inline map<ShowerParticlePtr,tNasonBranchingPtr> & particles();

  /**
   *  Access the set of branchings
   */
  inline set<NasonBranchingPtr> & branchings();
  
  /**
   * Access the incoming branchings
   */
  inline set<NasonBranchingPtr> & incoming();
private:

  /**
   *  The ShowerTree
   */
  ShowerTreePtr _tree;

  /**
   *  Map from the particles in the ShowerTree to the NasonBranchings
   */
  map<ShowerParticlePtr,tNasonBranchingPtr> _particles;

  /**
   *  The NasonBranchings in the hard process
   */
  set<NasonBranchingPtr> _branchings;

  /**
   *  The NasonBranchings which initiate the space-like showers
   */
  set<NasonBranchingPtr> _spacelike;
};

/**
 * The NasonBranching class is designed to contain the information needed for
 * an individual branching in the Nason approach
 */
class NasonBranching : public Base {

  /**
   *  The NasonTree is friend
   */
  friend class NasonTree;

public:

  /**
   * The default constructor
   * @param particle The particle which is branching
   * @param sudakov  The Sudakov form factor for the branching
   * @param parent   The parent for the branching
   * @param incoming Whether the particle is incoming or outgoing
   */
  inline NasonBranching(ShowerParticlePtr particle,
			SudakovPtr sudakov,
			tNasonBranchingPtr parent,bool incoming);

  /**
   * Add a child of the branching
   * @param child The child of the branching
   */
  inline void addChild(NasonBranchingPtr child);

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
  tNasonBranchingPtr _parent;

  /**
   *  The Sudakov form-factor
   */
  SudakovPtr _sudakov;

  /**
   * The children
   */
  vector<NasonBranchingPtr> _children;

  /**
   *  The beam particle
   */
  PPtr _beam;
};
}

#include "NasonTree.icc"

#endif /* HERWIG_NasonTree_H */
