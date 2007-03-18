// -*- C++ -*-
#ifndef HERWIG_NasonTree_H
#define HERWIG_NasonTree_H
//
// This is the declaration of the NasonTree class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
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
  NasonTree(vector<NasonBranchingPtr>);

private:

  /**
   *  The ShowerTree
   */
  ShowerTreePtr _tree;

  /**
   *  Map from the particles in the ShowerTree to the NasonBranchings
   */
  map<ShowerProgenitorPtr,NasonBranchingPtr> _particles;

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
   * @param children The branching products
   */
  inline NasonBranching(ShowerParticlePtr particle,bool incoming=false,
			vector<NasonBranchingPtr> children=vector<NasonBranchingPtr>());
  
private:

  /**
   *  The branching particle
   */
  ShowerParticlePtr _particle;

  /**
   *  The rescaled momentum
   */
  Lorentz5Momentum original;

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
   * The children
   */
  vector<NasonBranchingPtr> _children;
};
}

#include "NasonTree.icc"

#endif /* HERWIG_NasonTree_H */
