// -*- C++ -*-
#ifndef HERWIG_ShowerTree_H
#define HERWIG_ShowerTree_H

#include "ThePEG/Config/ThePEG.h"
#include "ShowerConfig.h"
#include "Herwig++/Shower/Kinematics/ShowerParticle.h"
#include "ShowerProgenitor.fh"
#include "ThePEG/EventRecord/Step.h"
#include <cassert>
#include "ShowerTree.fh"

namespace Herwig {
  
using namespace ThePEG;
 
/** \ingroup Shower
 * 
 *  The ShowerTree class stores the basic information needed for
 *  each hard interaction, either a scattering process or decay, which 
 *  needs to be showered.
 *
 */
class ShowerTree : public Base
{
public:

  /**
   *  Constructors and Destructors
   */
  //@{
  /**
   * Constructor for a scattering process
   * @param eh The event handler
   * @param out The outgoing particles
   * @param vars Pointer to the ShowerVariables object to provide access to some members
   * @param decay Map into which the trees for any unstable particles are inserted
   * @param ch Access to the event handler
   */
  ShowerTree(tEHPtr eh, const ParticleVector & out,
	     ShowerVarsPtr vars,
	     multimap<Energy,ShowerTreePtr> & decay);
  
  /**
   *  Constructor for a decay
   * @param in The decaying particle
   * @param vars Pointer to the ShowerVariables object to provide access to some members
   * @param decay Map into which the trees for any unstable particles are inserted
   * @param ch Access to the event handler
   */
  ShowerTree(PPtr in, ShowerVarsPtr vars,multimap<Energy,ShowerTreePtr> & decay,
	     tEHPtr ch);
  //@}

public:

  /**
   * Insert the tree into the event record
   * @param pstep The step into which the particles should be inserted
   * @param ISR Whether or not ISR is switched on
   * @param FSR Whether or not FSR is switched on
   */
  inline void fillEventRecord(StepPtr pstep,bool ISR,bool FSR);

  /**
   * Set the parent tree to this tree for trees which come from this one.
   * This needs to be run after the constructor.
   */
  void setParents();

  /**
   *  Perform the decay for a tree starting with an unstable particle
   *  @param decay The map of widths and ShowerTrees for the decays so that
   *  any unstable decay products can be added.
   * @param ch Access to the event handler
   */
  void decay(multimap<Energy,ShowerTreePtr> & decay,tEHPtr ch);

  /**
   * Access methods for the type of interaction
   */
  //@{
  /**
   *  Whether or not this is a scattering process
   */
  bool isHard() const;

  /**
   *  Whether or not this is a decay.
   */
  bool isDecay() const;
  //@}

  /**
   *  Flags relating to the application of the hard matrix element correction
   */
  //@{
  /**
   *  Was the hard matrix element correction applied
   */
  inline bool hardMatrixElementCorrection() const;

  /**
   *  Set whether or not the hard matrix element correction was applied
   */ 
  inline void hardMatrixElementCorrection(bool in);
  //@}

  /**
   *  Get the incoming shower particles
   */
  inline map<ShowerProgenitorPtr,tShowerParticlePtr> & incomingLines();

  /**
   *  Get the outgoing shower particles
   */
  inline map<ShowerProgenitorPtr,tShowerParticlePtr> & outgoingLines();

  /**
   *  Update the shower product for a final-state particle
   */
  void updateFinalStateShowerProduct(ShowerProgenitorPtr progenitor,
				     ShowerParticlePtr parent,
				     const ShowerParticleVector & children);

  /**
   *  Update the shower product for an initial-state particle
   */
  void updateInitialStateShowerProduct(ShowerProgenitorPtr progenitor,
				       ShowerParticlePtr newParent);

  /**
   *  Get the current final shower product for a final-state particle
   */
  inline tShowerParticlePtr getFinalStateShowerProduct(ShowerProgenitorPtr progenitor);

  /**
   * Add a final-state branching. This method removes the parent of the branching
   * from the list of particles at the end of the shower and inserts the children
   * @param parent The parent for the branching
   * @param children The outgoing particles in the branching
   */
  inline void addFinalStateBranching(ShowerParticlePtr parent,
				     const ShowerParticleVector & children);

  /**
   *  Add an initial-state branching. This method removes the oldParent of the
   *  branching and inserts the result of the backward evolution and the 
   *  outgoing particle into the relevant lists.
   * @param oldParent The particle being backward evolved
   * @param newParent The initial-state particle resulting from the backward evolution
   * @param otherChild The final-state particle produced in the evolution.
   */
  inline void addInitialStateBranching(ShowerParticlePtr oldParent,
				       ShowerParticlePtr newParent,
				       ShowerParticlePtr otherChild);

  /**
   *  Access and set the flag for whether this tree has been showered
   */
  //@{
  /**
   *  Access the flag
   */
  inline bool hasShowered() const;

  /**
   *  Set the flag
   */
  inline void hasShowered(bool);
  //@}

  /**
   *  Access the parent tree
   */
  inline ShowerTreePtr parent() const;

  /**
   *  Clear all the shower products so that the event can be reshowered
   * if the first attempt fail
   */
  void clear();

  /**
   *  Reset the particles resulting from the shower to those which started
   *  the shower
   */
  void resetShowerProducts();

protected:

  /**
   * Functions to add the shower to the event record.
   */
  //@{
  /**
   * Insert a hard process
   * @param pstep The step into which the particles should be inserted
   * @param ISR Whether or not ISR is switched on
   * @param FSR Whether or not FSR is switched on
   */
  void insertHard(StepPtr pstep,bool ISR,bool FSR);

  /**
   * Insert a decay process
   * @param pstep The step into which the particles should be inserted
   * @param ISR Whether or not ISR is switched on
   * @param FSR Whether or not FSR is switched on
   */
  void insertDecay(StepPtr pstep,bool ISR,bool FSR);

  /**
   * Recursively add the final-state shower from the particle to the event record.
   * @param particle The final-state particle
   * @param step The step
   */
  void addFinalStateShower(PPtr particle, StepPtr step);

  /**
   *  Add the initial-state shwoer from the particle to the step
   * @param particle The final-state particle
   * @param step The step
   * @param addchildren Add the children of the particle
   */
  void addInitialStateShower(PPtr particle, StepPtr step, bool addchildren=true);

  /**
   * Update the colour information of a particle prior to insertion into the
   * event record.
   * @param particle The particle for which the colour is updated.
   */
  void updateColour(PPtr particle);
  //@}

  /**
   * Isolate the colour of the process from the rest of the event.
   * Called in the constructor
   * @param The original particles
   * @param The colour isolated copies
   */
  void colourIsolate(const vector<PPtr> & original, const vector<PPtr> & copy);

private:
  
  /**
   *  The incoming ShowerParticles connected to the interaction
   *  as the index of a map with the particle the shower backward evolves
   *  them to as the value
   */
  map<ShowerProgenitorPtr,tShowerParticlePtr> _incomingLines;

  /**
   *  The outgoing ShowerParticles connected to the interaction
   *  as the index of a map with the particle the shower
   *  evolves them to as the value
   */
  map<ShowerProgenitorPtr,tShowerParticlePtr> _outgoingLines;

  /**
   *  The outgoing ShowerParticles at the end of the final-state shower
   */
  set<tShowerParticlePtr> _forward;

  /**
   *  The incoming ShowerParticles at the end of the initial-state shower
   */
  set<tShowerParticlePtr> _backward;

  /**
   *  Was the hard matrix element correction applied
   */
  bool _hardMECorrection;

  /**
   *  Was this a scattering process or a decay
   */
  bool _wasHard;

  /**
   *  Map of colour lines used to reset colours when inserted into the event
   */
  map<ColinePtr,ColinePtr> _colour;

  /**
   *  Map of particles in this Tree which are the initial particles in other
   *  trees
   */
  map<tShowerTreePtr,tShowerProgenitorPtr> _treelinks;

  /**
   *  The parent tree
   */
  tShowerTreePtr _parent;

  /**
   *  Pointer to the shower variables
   */
  ShowerVarsPtr _showerVariables;

  /**
   *  Has this tree showered
   */
  bool _hasShowered;
};
}

#include "ShowerTree.icc"

#endif
