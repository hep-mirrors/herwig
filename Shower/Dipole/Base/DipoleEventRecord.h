// -*- C++ -*-
//
// DipoleEventRecord.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleEventRecord_H
#define HERWIG_DipoleEventRecord_H
//
// This is the declaration of the DipoleEventRecord class.
//

#include "Herwig/Shower/ShowerEventRecord.h"
#include "Herwig/Shower/PerturbativeProcess.h"
#include "ThePEG/PDF/PDF.h"
#include "Dipole.h"
#include "DipoleChain.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Johannes Bellm
 *
 * \brief Generalized dipole splitting info to deal with subleading-N
 * splittings.
 */
class SubleadingSplittingInfo
  : public DipoleSplittingInfo {

public:

  /**
   * Default constructor
   */
  SubleadingSplittingInfo()
    : DipoleSplittingInfo() {}

  /**
   * Get the iterator of the emitter dipole chain
   */
  list<DipoleChain>::iterator emitterChain() const { return theEmitterChain; }

  /**
   * Get the iterator of the emitter dipole
   */
  list<Dipole>::iterator emitterDipole() const { return theEmitterDipole; }

  /**
   * Get the iterator of the spectator dipole chain
   */
  list<DipoleChain>::iterator spectatorChain() const { return theSpectatorChain; }

  /**
   * Get the iterator of the spectator dipole
   */
  list<Dipole>::iterator spectatorDipole() const { return theSpectatorDipole; }

  /**
   * Get the starting scale
   */
  Energy startScale() const { return theStartScale; }

  /**
   * Set the iterator of the emitter dipole chain
   */
  void emitterChain(list<DipoleChain>::iterator it) { theEmitterChain = it; }

  /**
   * Set the iterator of the emitter dipole
   */
  void emitterDipole(list<Dipole>::iterator it) { theEmitterDipole = it; }

  /**
   * Set the iterator of the spectator dipole chain
   */
  void spectatorChain(list<DipoleChain>::iterator it) { theSpectatorChain = it; }

  /**
   * Set the iterator of the spectator dipole
   */
  void spectatorDipole(list<Dipole>::iterator it) { theSpectatorDipole = it; }

  /**
   * Set the starting scale
   */
  void startScale(Energy s) { theStartScale = s; }

private:

  /**
   * Iterator of the emitter dipole chain
   */
  list<DipoleChain>::iterator theEmitterChain;

  /**
   * Iterator of the emitter dipole
   */
  list<Dipole>::iterator theEmitterDipole;

  /**
   * Iterator of the spectator dipole chain
   */
  list<DipoleChain>::iterator theSpectatorChain;

  /**
   * Iterator of the spectator dipole
   */
  list<Dipole>::iterator theSpectatorDipole;

  /**
   * The starting scale
   */
  Energy theStartScale;

};

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Stephen Webster
 *
 * \brief The DipoleEventRecord class is 
 * used internally by the dipole shower.
 */
class DipoleEventRecord : public ShowerEventRecord {

public:

  /**
   * The default constructor.
   */
  DipoleEventRecord() {}

  /**
   * The default destructor just cleans up.
   */
  ~DipoleEventRecord() { clear(); }

public:

  /**
   * Return any non-coloured outgoing particles in the
   * current subprocess.
   */
  PList& hard() { return theHard; }

  /**
   * Return any non-coloured outgoing particles in the
   * current subprocess.
   */
  const PList& hard() const { return theHard; }

  /**
   * Return the momentum of the hard system
   */
  const Lorentz5Momentum& pX() const { return thePX; }

  /**
   * Transform all intermediate, hard and outgoing
   * partciles using the given transformation.
   */
  void transform(const SpinOneLorentzRotation& rot);

public:

  /**
   * Return the dipole chains to be showered.
   */
  const list<DipoleChain>& chains() const { return theChains; }

  /**
   * Access the dipole chains to be showered.
   */
  list<DipoleChain>& chains() { return theChains; }

  /**
   * Return the dipole chains which ceased evolving.
   */
  const list<DipoleChain>& doneChains() const { return theDoneChains; }

  /**
   * Access the dipole chains which ceased evolving.
   */
  list<DipoleChain>& doneChains() { return theDoneChains; }

  /**
   * Return true, if there are chains to be
   * showered.
   */
  bool haveChain() const { return !theChains.empty(); }

  /**
   * Return the current dipole chain
   */
  DipoleChain& currentChain() { assert(haveChain()); return theChains.front(); }

  /**
   * Pop the current dipole chain
   */
  void popChain();

  /**
   * Remove the given chain.
   */
  void popChain(list<DipoleChain>::iterator);

  /**
   * Remove the given chains.
   */
  void popChains(const list<list<DipoleChain>::iterator>&);

  /**
   * Create a merged dipole index given two independent dipoles;
   * the first dipole is to provide the emitter.
   */
  DipoleIndex 
  mergeIndex(list<Dipole>::iterator firstDipole, const pair<bool,bool>& whichFirst,
	     list<Dipole>::iterator secondDipole, const pair<bool,bool>& whichSecond) const;

  /**
   * Create a SubleadingSplitingInfo given two independent dipoles;
   * the first dipole is to provide the emitter.
   */
  SubleadingSplittingInfo 
  mergeSplittingInfo(list<DipoleChain>::iterator firstChain, list<Dipole>::iterator firstDipole, 
		     const pair<bool,bool>& whichFirst,
		     list<DipoleChain>::iterator secondChain, list<Dipole>::iterator secondDipole, 
		     const pair<bool,bool>& whichSecond) const;

  /**
   * Return a list of all possible subleading-N emitting pairs
   */
  void getSubleadingSplittings(list<SubleadingSplittingInfo>&);

public:

  /**
   * Split the dipole pointed to by the given iterator.
   * Return references to the affected chains, and update
   * iterators pointing to the children in the returned
   * chains.
   */

  void split(list<Dipole>::iterator dip,
	     DipoleSplittingInfo& dsplit,
	     pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
	     DipoleChain*& firstChain, DipoleChain*& secondChain) {
    split(dip,theChains.begin(),dsplit,childIterators,firstChain,secondChain,false);
  }

  /**
   * Split the dipole pointed to by the given iterator
   * in the indicated chain, indicating a splitting with
   * a colour spectator.
   * Return references to the affected chains, and update
   * iterators pointing to the children in the returned
   * chains.
   */
  void split(list<Dipole>::iterator dip,
	     list<DipoleChain>::iterator ch,
	     DipoleSplittingInfo& dsplit,
	     pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
	     DipoleChain*& firstChain, DipoleChain*& secondChain,
	     bool colourSpectator = true);
  
  /**
   * As split, but not touching the acctual event record.
   */
  
  pair<PVector,PVector> tmpsplit(list<Dipole>::iterator dip,
             DipoleSplittingInfo& dsplit,
             pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
             DipoleChain*& firstChain, DipoleChain*& secondChain) {
    return tmpsplit(dip,theChains.begin(),dsplit,childIterators,firstChain,secondChain,false);
  }
  
  /**
   * As split, but not touching the acctual event record.
   */
  pair<PVector,PVector> tmpsplit(list<Dipole>::iterator dip,
             list<DipoleChain>::iterator ch,
             DipoleSplittingInfo& dsplit,
             pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
             DipoleChain*& firstChain, DipoleChain*& secondChain,
             bool colourSpectator = true);


  /**
   * Let the given dipole take the recoil of 
   * the indicated splitting.
   */
  void recoil(list<Dipole>::iterator dip,
	      list<DipoleChain>::iterator ch,
	      DipoleSplittingInfo& dsplit);

  /**
   * Peform a subleading-N splitting
   */
  void splitSubleading(SubleadingSplittingInfo& dsplit,
		       pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
		       DipoleChain*& firstChain, DipoleChain*& secondChain);

  /**
   * Update the particles upon insertion of the
   * given splitting.
   */
  void update(DipoleSplittingInfo& dsplit);
  
  /**
   * As update, but not touching the acctual event record.
   */
  pair<PVector,PVector> tmpupdate(DipoleSplittingInfo& dsplit);

  /**
   * Return the dipole(s) containing the incoming
   * partons after the evolution has ended. Put back
   * the chains containing these to the chains to be
   * showered.
   */
  list<pair<list<Dipole>::iterator,list<DipoleChain>::iterator> >
  inDipoles();

  /**
   * Fill the given step and return incoming partons.
   */
  tPPair fillEventRecord(StepPtr step, bool firstInteraction, bool realigned);

public:

  /**
   * Prepare the event record for the given
   * subprocess.
   */
  const map<PPtr,PPtr>& prepare(tSubProPtr subpro,
                                tStdXCombPtr xc,
				StepPtr step,
                                const pair<PDF,PDF>& pdf,
				tPPair beam,
				bool firstInteraction,
                                bool dipoles = true);
  /**
   * Prepare the event record for the given
   * subprocess.
   */
  void slimprepare(tSubProPtr subpro,
		   tStdXCombPtr xc,
		   const pair<PDF,PDF>& pdf,tPPair beam,
		   bool dipoles = true);

  /**
   * Clear the event record: Give up ownership
   * on any object involved in the evolution.
   */
  virtual void clear();

public:

  /**
   * Print event record at current state.
   */
  void debugLastEvent(ostream&) const;

public:

  /**
   *  Get the decays
   */
  map<PPtr,PerturbativeProcessPtr> & decays() {return theDecays;}

  /**
   * Used in DipoleEventRecord::prepare.
   * Add the outgoing particles from a perturbative 
   * process to the vector of original particles. 
   * Iterates through decay chains.
   **/
  void fillFromDecays(PerturbativeProcessPtr decayProc, vector<PPtr>& original);

  /**
   * Used in DipoleEventRecord::prepare.
   * Replace the particles in the given 
   * perturbative process with their copies from the 
   * map, theOriginals.
   * Iterates through decays chains.
   **/
  void separateDecay(PerturbativeProcessPtr decayProc);

  /**
   *  Decay the particle
   */
  Energy decay(PPtr incoming, bool& powhegEmission);

  /**
   * Prepare the event record for the showering of a decay.
   * Return false if the decay does not need to be showered.
   **/
  bool prepareDecay(PerturbativeProcessPtr decayProc);

  /**
   * Boost the momentum of the outgoing of the given 
   * perturbative process to the momentum of given particle.
   **/
  void updateDecayMom(PPtr decayParent, PerturbativeProcessPtr decayProc);

  /**
   * Iteratively update the momenta of all
   * particles in a decay chain, starting 
   * with the outgoing from the given parent
   **/
  void updateDecayChainMom(PPtr decayParent, PerturbativeProcessPtr decayProc);

  /**
   * Update theDecays following the decay and/or
   * showering of a decay particle.
   * With iteration switched on (true) this will
   * update theDecays with the entire decay chain.
   `**/
  void updateDecays(PerturbativeProcessPtr decayProc, bool iterate = true);

  /**
   *  Access current decay process
   */
  PerturbativeProcessPtr currentDecay() {return theCurrentDecay;}

  /**
   *  Set current decay process
   */
  void currentDecay(PerturbativeProcessPtr in) {theCurrentDecay=in;}


  // SW - Changed from protected to public so that functions can be used in DipoleShowerHandler
public:

  /**
   * Find the chains to be showered.
   * The decay bool avoids mixing up decaying particles in hard and decay processes
   */
  void findChains(const PList& ordered, const bool decay = false);

  /**
   * Sort the coloured partons into a colour ordered ensemble.
   */
  PList colourOrdered(PPair & in,PList & out);


private:

  struct getMomentum {
    const Lorentz5Momentum& operator() (PPtr particle) const {
      return particle->momentum();
    }
  };

  /**
   * The momentum of the hard system
   */
  Lorentz5Momentum thePX;

  /**
   * Any non-coloured outgoing particles in the
   * current subprocess.
   */
  PList theHard;

  /**
   * Map originals to copies.
   */
  map<PPtr,PPtr> theOriginals;

  /**
   * The dipole chains currently showered.
   */
  list<DipoleChain> theChains;

  /**
   * The dipole chains which ceased evolving.
   */
  list<DipoleChain> theDoneChains;

private:

  /**
   * Storage of the particles which need to be decayed
   */
  map<PPtr,PerturbativeProcessPtr> theDecays;

  /**
   *
   */
  PerturbativeProcessPtr theCurrentDecay;
};


}

#endif /* HERWIG_DipoleEventRecord_H */
