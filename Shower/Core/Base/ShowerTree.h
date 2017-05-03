// -*- C++ -*-
//
// ShowerTree.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerTree_H
#define HERWIG_ShowerTree_H

#include "ThePEG/Config/ThePEG.h"
#include "Herwig/Shower/ShowerHandler.fh"
#include "Herwig/Shower/PerturbativeProcess.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/ShowerEventRecord.h"
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "ShowerProgenitor.h"
#include "ThePEG/EventRecord/Step.h"
#include <cassert>
#include "ShowerTree.fh"

namespace Herwig {

/**
 *  Typedef for map of ShowerTrees for decays
 */
typedef multimap<Energy,ShowerTreePtr,std::greater<Energy> > ShowerDecayMap;
  
using namespace ThePEG;
 
/** \ingroup Shower
 * 
 *  The ShowerTree class stores the basic information needed for
 *  each hard interaction, either a scattering process or decay, which 
 *  needs to be showered.
 *
 */
class ShowerTree : public ShowerEventRecord {

  friend class ShowerHandler;

public:

  /**
   * Constructor from a perturbative process
   * @param process The perturbative process
   */
  ShowerTree(PerturbativeProcessPtr process);

  /**
   * Calculate the space-time displacement
   * @param particle The particle for which to calculate the displacement
   */
  static Lorentz5Distance spaceTimeDistance(tPPtr particle);

  /**
   *  Construct the trees from the hard process
   * @param hardTree The output ShowerTree for the hard process
   * @param decayTrees The output ShowerTrees for any decays.
   * @param hard The output ShowerTree for the hard process
   * @param decay The output ShowerTrees for any decays.
   */
  static void constructTrees(ShowerTreePtr & hardTree,
			     ShowerDecayMap & decayTrees,
			     PerturbativeProcessPtr hard,
			     DecayProcessMap decay);

public:

  /**
   * Insert the tree into the event record
   * @param pstep The step into which the particles should be inserted
   * @param ISR Whether or not ISR is switched on
   * @param FSR Whether or not FSR is switched on
   */
  void fillEventRecord(StepPtr pstep,bool ISR,bool FSR) {
    if(_wasHard) 
      insertHard (pstep,ISR,FSR);
    else         
      insertDecay(pstep,ISR,FSR);
  }
  

  /**
   * Set the parent tree to this tree for trees which come from this one.
   * This needs to be run after the constructor.
   */
  void setParents();

  /**
   * Access methods for the type of interaction
   */
  //@{
  /**
   *  Whether or not this is a scattering process
   */
  bool isHard() const { return _wasHard; } 


  /**
   *  Whether or not this is a decay.
   */
  bool isDecay() const { return !_wasHard; }
  //@}

  /**
   *  Flags relating to the application of the hard matrix element correction
   */
  //@{
  /**
   *  Was the hard matrix element correction applied
   */
  bool hardMatrixElementCorrection() const { return _hardMECorrection; }

  /**
   *  Set whether or not the hard matrix element correction was applied
   */ 
  void hardMatrixElementCorrection(bool in) { _hardMECorrection=in; }
  //@}

  /**
   *  Get the incoming shower particles
   */
  map<ShowerProgenitorPtr,ShowerParticlePtr> & incomingLines() {
    return _incomingLines; 
  }

  /**
   *  Get the outgoing shower particles
   */
  map<ShowerProgenitorPtr,tShowerParticlePtr> & outgoingLines() {
    return _outgoingLines; 
  }

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
  tShowerParticlePtr getFinalStateShowerProduct(ShowerProgenitorPtr progenitor) {
    return _outgoingLines.find(progenitor)==_outgoingLines.end()
      ? tShowerParticlePtr() : _outgoingLines[progenitor];
  }

  /**
   * Add a final-state branching. This method removes the parent of the branching
   * from the list of particles at the end of the shower and inserts the children
   * @param parent The parent for the branching
   * @param children The outgoing particles in the branching
   */
  void addFinalStateBranching(ShowerParticlePtr parent,
			      const ShowerParticleVector & children);

  /**
   *  Add an initial-state branching. This method removes the oldParent of the
   *  branching and inserts the result of the backward evolution and the 
   *  outgoing particle into the relevant lists.
   * @param oldParent The particle being backward evolved
   * @param newParent The initial-state particle resulting from the backward evolution
   * @param otherChild The final-state particle produced in the evolution.
   */
  void addInitialStateBranching(ShowerParticlePtr oldParent,
				ShowerParticlePtr newParent,
				ShowerParticlePtr otherChild);

  // /**
  //  *  Member called at the end of the shower of a tree to perform a number of
  //  *  updates.
  //  *  @param decay The map of widths and ShowerTrees for the decays so that
  //  *  any unstable decay products can be added.
  //  */
  void updateAfterShower(ShowerDecayMap & decay);

  /**
   *  Access and set the flag for whether this tree has been showered
   */
  //@{
  /**
   *  Access the flag
   */
  bool hasShowered() const { return _hasShowered; }

  /**
   *  Set the flag
   */
  void hasShowered(bool in) { _hasShowered=in; }
  //@}

  /**
   *  Access the parent tree
   */
  ShowerTreePtr parent() const { return _parent; }

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

  /**
   *  Set maximum Emission scales
   */
  void setVetoes(const map<ShowerInteraction,Energy> & pTs,
		 unsigned int type);

  /**
   *  Extract the progenitors for the reconstruction
   */
  vector<ShowerProgenitorPtr> extractProgenitors();

  /**
   *  Access to the outgoing particles
   */
  const set<tShowerParticlePtr> & forwardParticles() const { return _forward; }

  /**
   *  Map of particles in this Tree which are the initial particles in other
   *  trees
   */
  const map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> > &
  treelinks()  const {return _treelinks;}

  /**
   *  Update the link between shower particle and tree
   */
  void updateLink(tShowerTreePtr tree,
		  pair<tShowerProgenitorPtr,tShowerParticlePtr> in) {
    _treelinks[tree] = in;
  }

  /**
   *  Transform the tree
   */
  void transform(const LorentzRotation & rot, bool applyNow);

  /**
   *  Apply any postphoned transformations
   */
  void applyTransforms();

  /**
   *   Clear any postphoned transformations
   */ 
  void clearTransforms();

  /**
   *  Transform which needs to be applied
   */
  const LorentzRotation & transform() {return _transforms;}

  /**
   *  Get all the progenitors
   */
  vector<ShowerParticlePtr> extractProgenitorParticles();

  /**
   *    Check the momentum conservation in the tree
   */
  void checkMomenta();

  /**
   *  Update tree after the parent has been decayed.
   */
  void update(PerturbativeProcessPtr newProcess);

  /**
   *  The perturbative process
   */
  RealEmissionProcessPtr perturbativeProcess();

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
   * @param hadron The incoming hadron
   * @param step The step
   * @param addchildren Add the children of the particle
   */
  void addInitialStateShower(PPtr particle, PPtr hadron,
			     StepPtr step, bool addchildren=true);
  //@}

  /**
   *  After the creatation of a ShowerParticle make sure it is properly attached 
   *  to its ColourLine
   * @param part The particle
   */
  void fixColour(tShowerParticlePtr part);

private:

  /**
   * Incoming partons for the hard process
   */
  PPair _incoming;
  
  /**
   *  The incoming ShowerParticles connected to the interaction
   *  as the index of a map with the particle the shower backward evolves
   *  them to as the value
   */
  map<ShowerProgenitorPtr,ShowerParticlePtr> _incomingLines;

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
   *  Map of particles in this Tree which are the initial particles in other
   *  trees
   */
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> > _treelinks;

  /**
   *  The parent tree
   */
  tShowerTreePtr _parent;

  /**
   *  Has this tree showered
   */
  bool _hasShowered;

  /**
   *  The transforms which still need to be applied
   */
  LorentzRotation _transforms;

private:

  /**
   *  Whether or not to include space-time distances
   */
  static bool _spaceTime;

  /**
   *  Minimum virtuality for the space-time model
   */
  static Energy2 _vmin2;

};
}

#endif
