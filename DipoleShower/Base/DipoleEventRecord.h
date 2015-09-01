// -*- C++ -*-
//
// DipoleEventRecord.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleEventRecord_H
#define HERWIG_DipoleEventRecord_H
//
// This is the declaration of the DipoleEventRecord class.
//

#include "ThePEG/PDF/PDF.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Dipole.h"
#include "DipoleChain.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
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
 * \author Simon Platzer
 *
 * \brief The DipoleEventRecord class is 
 * used internally by the dipole shower.
 */
class DipoleEventRecord {

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
   * Return the incoming partons at the current 
   * stage of the evolution.
   */
  PPair& incoming() { return theIncoming; }

  /**
   * Return the incoming partons at the current 
   * stage of the evolution.
   */
  const PPair& incoming() const { return theIncoming; }

  /**
   * Return the outgoing partons at the current
   * stage of the evolution.
   */
  PList& outgoing() { return theOutgoing; }

  /**
   * Return the outgoing partons at the current
   * stage of the evolution.
   */
  const PList& outgoing() const { return theOutgoing; }

  /**
   * Return the intermediate particles at the current
   * stage of the evolution.
   */
  PList& intermediates() { return theIntermediates; }

  /**
   * Return the intermediate particles at the current
   * stage of the evolution.
   */
  const PList& intermediates() const { return theIntermediates; }

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
   * Return the subprocess currently showered
   */
  tSubProPtr subProcess() const { return theSubProcess; }

  /**
   * Return the XComb describing the hard process.
   */
  tStdXCombPtr xcombPtr() const { return theXComb; }

  /**
   * Return the XComb describing the hard process.
   */
  const StandardXComb& xcomb() const { return *theXComb; }

  /**
   * Return the momentum fractions.
   */
  const pair<double,double>& fractions() const { return theFractions; }

  /**
   * Return the PDFs
   */
  const pair<PDF,PDF>& pdfs() const { return thePDFs; }

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
				const pair<PDF,PDF>& pdf,
				bool dipoles = true);

  /**
   * Clear the event record: Give up ownership
   * on any object involved in the evolution.
   */
  void clear();

public:

  /**
   * Print event record at current state.
   */
  void debugLastEvent(ostream&) const;

protected:

  /**
   * Sort the coloured partons into a colour ordered ensemble.
   */
  PList colourOrdered();

  /**
   * Find the chains to be showered.
   */
  void findChains(const PList& ordered);

  /**
   * Add all particles to the relevant sets
   */
  void getAll(const ParticleVector& childs,
	      set<PPtr>& hardSet,
	      set<PPtr>& outgoingSet);

  /**
   * Isolate the colour of the process from the rest of the event.
   * Called in the constructor
   */
  void colourIsolate(const vector<PPtr> & original, const vector<PPtr> & copy);

  /**
   * Update the colour information of a particle prior to insertion into the
   * event record.
   */
  void updateColour(PPtr particle);

private:

  struct getMomentum {
    const Lorentz5Momentum& operator() (PPtr particle) const {
      return particle->momentum();
    }
  };

  /**
   * The subprocess currently showered.
   */
  SubProPtr theSubProcess;

  /**
   * Pointer to the XComb which generated the hard process.
   */
  StdXCombPtr theXComb;

  /**
   * The PDFs to be considered.
   */
  pair<PDF,PDF> thePDFs;

  /**
   * The momentum of the hard system
   */
  Lorentz5Momentum thePX;

  /**
   * Momentum fractions of the incoming partons.
   */
  pair<double,double> theFractions;

  /**
   * The incoming partons at the current
   * stage of the evolution.
   */
  PPair theIncoming;

  /**
   * The outgoing partons at the current
   * stage of the evolution.
   */
  PList theOutgoing;

  /**
   * The intermediate particles at the current
   * stage of the evolution.
   */
  PList theIntermediates;

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
   * Map colour lines from copies to originals.
   */
  map<ColinePtr,ColinePtr> theColourLines;

  /**
   * The dipole chains currently showered.
   */
  list<DipoleChain> theChains;

  /**
   * The dipole chains which ceased evolving.
   */
  list<DipoleChain> theDoneChains;

};


}

#endif /* HERWIG_DipoleEventRecord_H */
