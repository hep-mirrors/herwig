  // -*- C++ -*-
  //
  // DipoleShowerHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the DipoleShowerHandler class.
  //

#include <config.h>
#include "DipoleShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

  // include theses to have complete types
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "Herwig/PDF/HwRemDecayer.h"

#include "Herwig/Shower/Dipole/Utility/DipolePartonSplitter.h"
#include "Herwig/MatrixElement/Matchbox/Base/MergerBase.h"

#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include <queue>

using namespace Herwig;

bool DipoleShowerHandler::firstWarn = true;

DipoleShowerHandler::DipoleShowerHandler() :
  ShowerHandler(), chainOrderVetoScales(true),
  nEmissions(0), discardNoEmissions(false), firstMCatNLOEmission(false),
  thePowhegDecayEmission(true),
  //theAnalyseSpinCorrelations(false),
  realignmentScheme(0),
  doSubleadingNc(false),subleadingNcEmissionsLimit(0),
  densityOperatorEvolution(0),densityOperatorCutoff(1.0*GeV2),
  doPartialUnweightingAtEmission(false),
  doPartialUnweighting(false),referenceWeight(0.1),
  cmecReweightFactor(1.0),negCMECScaling(1.0),
  verbosity(0), printEvent(0), nTries(0),
  didRadiate(false), didRealign(false),
  theRenormalizationScaleFreeze(1.*GeV),
  theFactorizationScaleFreeze(2.*GeV), theDoCompensate(false),
  theFreezeGrid(500000), theDetuning(1.0),
  maxPt(ZERO), muPt(ZERO),
  theInputColouredOffShellInShower(),
  theZBoundaries(1) {}

DipoleShowerHandler::~DipoleShowerHandler() {}

IBPtr DipoleShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleShowerHandler::fullclone() const {
  return new_ptr(*this);
}



void  DipoleShowerHandler::cascade(tPVector ) {
  throw Exception()
  << "DipoleShowerHandler: Dipoleshower not implemented as second shower."
  << "Check your setup or contact Herwig authors."
  << Exception::runerror;
}


tPPair DipoleShowerHandler::cascade(tSubProPtr sub, XCombPtr,
                                    Energy optHardPt, Energy optCutoff) {
  
  useMe();
  
  prepareCascade(sub);
  resetWeights();
  
  if ( !doFSR() && ! doISR() )
    return sub->incoming();

  eventRecord().setSubleadingNc(doSubleadingNc,
				subleadingNcEmissionsLimit);
  eventRecord().clear();
  eventRecord().prepare(sub,dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr()),newStep(),pdfs(),
			ShowerHandler::currentHandler()->generator()->currentEvent()->incoming(),
			firstInteraction(), offShellPartons(),
                        !doSubleadingNc);
  if ( doSubleadingNc ) {
    if ( !theSplittingReweight ) {
      throw Exception() << "No splitting reweight was found. "
			<< "A ColourMatrixElementCorrection "
			<< "splitting reweight is required "
			<< "for the subleading colour shower."
			<< Exception::runerror; 
    }
    //Set the evolution scheme for the density operator
    eventRecord().setDensityOperatorEvolution( densityOperatorEvolution, densityOperatorCutoff );
    //Set the CMEC reweight factor
    theSplittingReweight->reweightFactor(cmecReweightFactor);
    theSplittingReweight->negativeScaling(negCMECScaling);
    theSplittingReweight->updateCurrentHandler();
  }

  // SW: Removed simple test on doFSR and doISR and moved
  // here to account for the case of a hard event involving
  // no coloured particles but with unstable outgoing particles
  if ( !doFSR() && ! doISR() && eventRecord().decays().empty() )
  return sub->incoming();
  if ( !doISR() &&
       eventRecord().outgoing().empty() &&
       eventRecord().decays().empty() )
    return sub->incoming();
  if ( !doFSR() &&
       !eventRecord().incoming().first->coloured() &&
      !eventRecord().incoming().second->coloured() &&
       eventRecord().decays().empty() )
  return sub->incoming();
  
  nTries = 0;
  
  // Clear the vertex record for spin correlations
  if ( spinCorrelations() ) //|| theAnalyseSpinCorrelations )
    vertexRecord().clear();
  
  while ( true ) {
    
    try {
      
      didRadiate = false;
      didRealign = false;
      
      if ( eventRecord().truncatedShower() ) {
        throw Exception() << "Inconsistent hard emission set-up in DipoleShowerHandler::cascade. "
        << "No truncated shower needed with DipoleShowerHandler.  Add "
        << "'set MEMatching:TruncatedShower No' to input file."
        << Exception::runerror;
      }
      
      hardScales(lastXCombPtr()->lastShowerScale());
      
      if ( verbosity > 1 ) {
        generator()->log() << "DipoleShowerHandler starting off:\n";
        eventRecord().debugLastEvent(generator()->log());
        generator()->log() << flush;
      }
      
      unsigned int nEmitted = 0;
      
      if ( firstMCatNLOEmission ) {
        
        if ( !eventRecord().isMCatNLOHEvent() )
        nEmissions = 1;
        else
        nEmissions = 0;
        
      }
      
      if ( !firstMCatNLOEmission ) {
        
        doCascade(nEmitted,optHardPt,optCutoff);
        
        if ( discardNoEmissions ) {
          if ( !didRadiate )
          throw Veto();
          if ( nEmissions )
          if ( nEmissions < nEmitted )
          throw Veto();
        }
        
      } else {
        
        if ( nEmissions == 1 )
        doCascade(nEmitted,optHardPt,optCutoff);
        
      }
      
      if ( intrinsicPtGenerator ) {
        if ( eventRecord().incoming().first->coloured() &&
            eventRecord().incoming().second->coloured() ) {
          LorentzRotation rot =
          intrinsicPtGenerator->kick(eventRecord().incoming(),
                                     eventRecord().intermediates());
          eventRecord().transform(rot);
        }
      }
      
      didRealign = realign();
    
      constituentReshuffle();

      // backup subleading switch if decays fail
      bool doneSubleadingNc = doSubleadingNc;

      // subleading N can't handle decays
      doSubleadingNc = false;

      try {
      
        // Decay and shower any particles that require decaying
      while ( !eventRecord().decays().empty() ) {
	
        map<PPtr,PerturbativeProcessPtr>::const_iterator decayIt = eventRecord().decays().begin();
        if ( eventRecord().nextDecay() ) {
          decayIt = eventRecord().decays().find(eventRecord().nextDecay() );
        }
        else {
          // find the decay to do, one with greatest width and parent showered
        while(find(eventRecord().outgoing().begin(),eventRecord().outgoing().end(),decayIt->first)==
              eventRecord().outgoing().end() &&
              find(eventRecord().hard().begin(),eventRecord().hard().end(),decayIt->first)==
              eventRecord().hard().end()) ++decayIt;
        }
	
        assert(decayIt!=eventRecord().decays().end());
        PPtr incoming = decayIt->first;
        eventRecord().currentDecay(decayIt->second);
        
          // Use this to record if an emission actually happens
        bool powhegEmission = !( nEmissions && nEmitted==nEmissions) ? thePowhegDecayEmission : false;
        
          // Decay the particle / sort out its pert proc
        Energy showerScale = eventRecord().decay(incoming, powhegEmission);
        
          // Following the decay, the bool powheg emission is updated
          // to indicate whether or not an emission occurred
        if ( powhegEmission )
        nEmitted += 1;
        
          // Check that there is only one particle incoming to the decay
        assert(eventRecord().currentDecay()->incoming().size()==1);
        
          // Prepare the event record for the showering of the decay
        bool needToShower = eventRecord().prepareDecay(eventRecord().currentDecay(),
						       offShellPartons());
        
          // Only need to shower if we have coloured outgoing particles
        if ( needToShower ) {
          
            // The decays currently considered produce a maximum of 2 chains (with powheg emission)
            // so all dipole should have the same scale as returned by the decay function.
          assert( eventRecord().chains().size() <= 2 );
          for ( auto  & ch : eventRecord().chains()) {
            for ( auto  & dip : ch.dipoles()) {
              assert ( showerScale > ZERO );
              dip.leftScale( showerScale );
              dip.rightScale( showerScale );
            }
          }
      
          // Prepare vertex record for spin correlations in decay shower
          if ( spinCorrelations() )
            vertexRecord().prepareParticleDecay(incoming);
      
            // Perform the cascade
          doCascade(nEmitted,optHardPt,optCutoff,true);
          
          if ( spinCorrelations() )
            vertexRecord().updateParticleDecay();
      
            // Do the constituent mass shell reshuffling
          decayConstituentReshuffle(eventRecord().currentDecay());
          
        }
        
          // Update the decays, adding any decays and updating momenta
        eventRecord().updateDecays(eventRecord().currentDecay());
        
	  eventRecord().decays().erase(decayIt);

	}

      } catch(...) {

	// reset flag
	doSubleadingNc = doneSubleadingNc;
	throw;

      }

      doSubleadingNc = doneSubleadingNc;
      
      break;
      
    } catch (RedoShower&) {
      
      resetWeights();
      
      if ( ++nTries > maxtry() )
      throw ShowerTriesVeto(maxtry());
      
      eventRecord().clear();
      // reset the spininfos
      if(sub->incoming().first->spinInfo())
	sub->incoming().first->spinInfo()->reset();
      if(sub->incoming().second->spinInfo())
	sub->incoming().second->spinInfo()->reset();
      for(PPtr out : sub->outgoing()) {
	if(out->spinInfo()) out->spinInfo()->reset();
      }
      // prepare for the new shower
      eventRecord().prepare(sub, dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr()), newStep(), pdfs(),
                            ShowerHandler::currentHandler()->generator()->currentEvent()->incoming(),
			    firstInteraction(), offShellPartons(),
                            !doSubleadingNc);
      if ( doSubleadingNc ) {
	theSplittingReweight->updateCurrentHandler();
      }
      
      continue;
      
    } catch (...) {
      throw;
    }
    
  }
  
  tPPair incoming=eventRecord().fillEventRecord(newStep(),firstInteraction(),didRealign);
  setDidRunCascade(true);
  return incoming;
}



  // Reshuffle the outgoing partons from the hard process onto their constituent mass shells
void DipoleShowerHandler::constituentReshuffle() {
  
  if ( constituentReshuffler &&  ShowerHandler::currentHandler()->retConstituentMasses() ) {
    if ( eventRecord().decays().empty() ) {
      constituentReshuffler->reshuffle(eventRecord().outgoing(),
                                       eventRecord().incoming(),
                                       eventRecord().intermediates());
      return;
    }
    
    else {
      PList decaying;
      for(auto const & dec : eventRecord().decays())
      decaying.push_back(dec.first);
      
      constituentReshuffler->hardProcDecayReshuffle( decaying,
                                                    eventRecord().outgoing(),
                                                    eventRecord().hard(),
                                                    eventRecord().incoming(),
                                                    eventRecord().intermediates());
    }
  }
  
    // After reshuffling the hard process, the decays need to be updated
    // as this is not done in reshuffle
  vector<pair<PPtr,PerturbativeProcessPtr> > decays;
  for(auto const & dec : eventRecord().decays() )
  decays.push_back({dec.first,dec.second});
  
  
  for(auto const & dec : decays) {
    
    PPtr unstable = dec.first;
    PList::iterator pos = find(eventRecord().intermediates().begin(),
                               eventRecord().intermediates().end(),
                               dec.first);
    
      // Update the PPtr in theDecays
    if(pos!=eventRecord().intermediates().end()) {
      unstable = *pos;
      while(!unstable->children().empty()) {
        unstable = unstable->children()[0];
      }
      eventRecord().decays().erase(dec.first);
      eventRecord().decays()[unstable] = dec.second;
      
        // Update the momenta of any other particles in the decay chain
        // (for externally provided events)
      if ( !(eventRecord().decays()[unstable]->outgoing().empty()) )
      eventRecord().updateDecayChainMom( unstable , eventRecord().decays()[unstable]);
    }
    
    else {
      if ( !(eventRecord().decays()[unstable]->outgoing().empty()) ) {
          // Update the momenta of any other particles in the decay chain
          // (for externally provided events)
          // Note this needs to be done for all decaying particles in the
          // outgoing/hard regardless of whether that particle radiated
          // or was involved in the reshuffling, this is due to the
          // transformation performed for IILightKinematics.
        if ( (find(eventRecord().outgoing().begin(),
                   eventRecord().outgoing().end(), unstable) != eventRecord().outgoing().end())
            || (find(eventRecord().hard().begin(),
                     eventRecord().hard().end(), unstable) != eventRecord().hard().end()) )
        eventRecord().updateDecayChainMom( unstable , eventRecord().decays()[unstable]);
        
      }
    }
  }
  
  eventRecord().currentDecay(PerturbativeProcessPtr());
}


  // Reshuffle outgoing partons from a decay process onto their constituent mass shells
void DipoleShowerHandler::decayConstituentReshuffle(PerturbativeProcessPtr decayProc) {
  
  
  
  if ( Debug::level  > 2  ){
  
      // Test this function by comparing the
      // invariant mass of the outgoing decay
      // systems before and after reshuffling
    Lorentz5Momentum testOutMomBefore (ZERO,ZERO,ZERO,ZERO);
    Energy testInvMassBefore = ZERO;
    
    for ( auto const & testDecayOutItBefore : decayProc->outgoing() ) {
      testOutMomBefore += testDecayOutItBefore.first->momentum();
    }
    
    testInvMassBefore = testOutMomBefore.m();
    
    
      // decayReshuffle updates both the event record and the decay perturbative process
    if ( constituentReshuffler && ShowerHandler::currentHandler()->retConstituentMasses()) {
      constituentReshuffler->decayReshuffle(decayProc,
                                            eventRecord().outgoing(),
                                            eventRecord().hard(),
                                            eventRecord().intermediates());
    }
    
    Lorentz5Momentum testOutMomAfter (ZERO,ZERO,ZERO,ZERO);
    Energy testInvMassAfter = ZERO;
    
    for ( auto const &  testDecayOutItAfter : decayProc->outgoing() ) {
      
      testOutMomAfter += testDecayOutItAfter.first->momentum();
      
    }
    
    testInvMassAfter = testOutMomAfter.m();
    
#ifndef NDEBUG
    Energy incomingMass = decayProc->incoming()[0].first->momentum().m();
#endif
    assert( abs(testInvMassBefore-incomingMass)/GeV < 1e-5 );
    assert( abs(testInvMassBefore-testInvMassAfter)/GeV < 1e-5);
    
    
  }else{
      // decayReshuffle updates both the event record and the decay perturbative process
    if ( constituentReshuffler && ShowerHandler::currentHandler()->retConstituentMasses() ) {
      constituentReshuffler->decayReshuffle(decayProc,
                                            eventRecord().outgoing(),
                                            eventRecord().hard(),
                                            eventRecord().intermediates());
    }
    return;
  }
  
}

  // Sets the scale of each particle in the dipole chains by finding the smallest
  //of several upper bound energy scales: the CMEnergy of the event,
  //the transverse mass of outgoing particles, the hardScale (maxPT or maxQ)
  //calculated for each dipole (in both configurations) and the veto scale for each particle
void DipoleShowerHandler::hardScales(Energy2 muf)  {
  
    // Initalise maximum pt as max CMEnergy of the event
  maxPt = generator()->maximumCMEnergy();
  
  if ( restrictPhasespace() ) {
      // First interaction == hard collision (i.e. not a MPI collision)
    if ( !hardScaleIsMuF() || !firstInteraction() ) {
      if ( !eventRecord().outgoing().empty() ) {
        for ( auto const &  p : eventRecord().outgoing() )
	       maxPt = min(maxPt,p->momentum().mt());
      }
        //Look at any non-coloured outgoing particles in the current subprocess
      else {
        assert(!eventRecord().hard().empty());
        Lorentz5Momentum phard(ZERO,ZERO,ZERO,ZERO);
        for ( auto const & p : eventRecord().hard())
        phard += p->momentum();
        Energy mhard = phard.m();
        maxPt = mhard;
      }
      maxPt *= hardScaleFactor();
    }
    else {
      maxPt = hardScaleFactor()*sqrt(muf);
    }
    muPt = maxPt;
  } else {
    muPt = hardScaleFactor()*sqrt(muf);
  }

  if ( doSubleadingNc ) {
    return;
  }

 for ( auto  & ch : eventRecord().chains()) {
    
      // Note that minVetoScale is a value for each DipoleChain, not each dipole
      // It will contain the minimum veto scale from all of the dipoles in the chain
    Energy minVetoScale = -1.*GeV;
    
    for ( auto  & dip : ch.dipoles()) {
      
        // max scale per config
      Energy maxFirst = ZERO;
      Energy maxSecond = ZERO;
      
        // Loop over the kernels for the given dipole.
        // For each dipole configuration, calculate ptMax (or QMax if virtuality ordering)
        // for each kernel and find the maximum
      for ( auto const & k : kernels) {
        
        pair<bool,bool> conf = {true,false};
        
        if ( k->canHandle(dip.index(conf)) ) {
            // Look in DipoleChainOrdering for this
          Energy scale =
          evolutionOrdering()->hardScale(dip.emitter(conf),dip.spectator(conf),
                                         dip.emitterX(conf),dip.spectatorX(conf),
                                         *k,dip.index(conf));
          maxFirst = max(maxFirst,scale);
        }
        
        conf = {false,true};
        
        if ( k->canHandle(dip.index(conf)) ) {
          Energy scale =
          evolutionOrdering()->hardScale(dip.emitter(conf),dip.spectator(conf),
                                         dip.emitterX(conf),dip.spectatorX(conf),
                                         *k,dip.index(conf));
          maxSecond = max(maxSecond,scale);
        }
        
      }
      
        // Find the maximum value from comparing the maxScale found from maxPt and the vetoScale of the particle
      if ( dip.leftParticle()->vetoScale() >= ZERO ) {
        maxFirst = min(maxFirst,sqrt(dip.leftParticle()->vetoScale()));
        
          // minVetoScale is a value for each DipoleChain, not each dipole
          // It contains the minimum veto scale for all the dipoles in the entire DipoleChain
        if ( minVetoScale >= ZERO )
        minVetoScale = min(minVetoScale,sqrt(dip.leftParticle()->vetoScale()));
        else
        minVetoScale = sqrt(dip.leftParticle()->vetoScale());
      }
      
      if ( dip.rightParticle()->vetoScale() >= ZERO ) {
        maxSecond = min(maxSecond,sqrt(dip.rightParticle()->vetoScale()));
        if ( minVetoScale >= ZERO )
        minVetoScale = min(minVetoScale,sqrt(dip.rightParticle()->vetoScale()));
        else
        minVetoScale = sqrt(dip.rightParticle()->vetoScale());
      }
      
        // Set the emitterScale for both members of each dipole
      maxFirst = min(maxPt,maxFirst);
      dip.emitterScale({true,false},maxFirst);
      
      maxSecond = min(maxPt,maxSecond);
      dip.emitterScale({false,true},maxSecond);
      
    }
    
      // if the smallest veto scale (i.e. from all of the dipoles)
      // is smaller than the scale calculated for a particular
      // particle in a particular dipole,
      // replace the scale with the veto scale
    if ( !evolutionOrdering()->independentDipoles() &&
        chainOrderVetoScales &&
        minVetoScale >= ZERO ) {
      for ( auto  & dip : ch.dipoles() ) {
        dip.leftScale(min(dip.leftScale(),minVetoScale));
        dip.rightScale(min(dip.rightScale(),minVetoScale));
      }
    }
    
  }
  
}

void DipoleShowerHandler::hardScalesSubleading(list<DipoleSplittingInfo> candidates,
					       Energy hardPt) {

  maxPt = hardPt;//generator()->maximumCMEnergy();

  // Note that minVetoScale is a value for each competing dipole (i.e. all dipoles
  // for the subleading shower.
  // It will contain the minimum veto scale from all of the dipoles
  Energy minVetoScale = -1.*GeV;

  for ( list<DipoleSplittingInfo>::iterator cand = candidates.begin();
	cand != candidates.end(); ++cand ) {

      // max scale
      Energy maxScale = ZERO;

      // Loop over kernels
      for ( vector<Ptr<DipoleSplittingKernel>::ptr>::iterator k =
	      kernels.begin(); k != kernels.end(); ++k ) {
	
	if ( (**k).canHandle(cand->index()) ) {
	  Energy scale =
	    evolutionOrdering()->hardScale(cand->emitter(),cand->spectator(),
					 cand->emitterX(),cand->spectatorX(),
					 **k,cand->index());
	  maxScale = max(maxScale,scale);
	}

      }

      if ( cand->emitter()->vetoScale() >= ZERO ) {
	maxScale = min(maxScale,sqrt(cand->emitter()->vetoScale()));
	if ( minVetoScale >= ZERO )
	  minVetoScale = min(minVetoScale,sqrt(cand->emitter()->vetoScale()));
	else
	  minVetoScale = sqrt(cand->emitter()->vetoScale());
      }

      maxScale = min(maxPt,maxScale);
      cand->scale(maxScale);

    }

    if ( !evolutionOrdering()->independentDipoles() &&
	 chainOrderVetoScales &&
	 minVetoScale >= ZERO ) {
      for ( list<DipoleSplittingInfo>::iterator cand = candidates.begin();
	    cand != candidates.end(); ++cand ) {
	cand->scale(min(cand->scale(),minVetoScale));
      }
    }

}



void DipoleShowerHandler::addCandidates(PPair particles,
					list<DipoleSplittingInfo>& clist) const {

  DipoleSplittingInfo candidate;
  Energy2 scale = ZERO;
  pair<bool,bool> is(particles.first == eventRecord().incoming().first,
		     particles.second == eventRecord().incoming().second);
  if ( (is.first && !is.second) ||
       (!is.first && is.second) ) {
    scale = -(particles.first->momentum() - particles.second->momentum()).m2();
  } else {
    scale = (particles.first->momentum() + particles.second->momentum()).m2();
  }

  DipoleIndex index(particles.first->dataPtr(),particles.second->dataPtr(),
		    is.first ? eventRecord().pdfs().first : PDF(),
		    is.second ? eventRecord().pdfs().second : PDF());

  candidate.scale(sqrt(scale));

  candidate.index(index);
  candidate.configuration(make_pair(true,false));
  candidate.emitter(particles.first);
  candidate.emitterX(is.first ? eventRecord().fractions().first : 1.0);
  candidate.spectator(particles.second);
  candidate.spectatorX(is.second ? eventRecord().fractions().second : 1.0);

  clist.push_back(candidate);

  index.swap();

  candidate.index(index);
  candidate.configuration(make_pair(false,true));
  candidate.emitter(particles.second);
  candidate.emitterX(is.second ? eventRecord().fractions().second : 1.0);
  candidate.spectator(particles.first);
  candidate.spectatorX(is.first ? eventRecord().fractions().first : 1.0);

  clist.push_back(candidate);

}

void DipoleShowerHandler::getCandidates(list<DipoleSplittingInfo>& clist) const {

  clist.clear();

  for ( PList::const_iterator i = eventRecord().outgoing().begin(); 
	i != eventRecord().outgoing().end(); ++i ) {
    PList::const_iterator j = i; ++j;
    for ( ; j != eventRecord().outgoing().end(); ++j ) {
      addCandidates(make_pair(*i,*j),clist);
    }
    // Changed order of *i and inc().first
    if ( eventRecord().incoming().first->coloured() )
      addCandidates(make_pair(eventRecord().incoming().first,*i),clist);
    if ( eventRecord().incoming().second->coloured() )
      addCandidates(make_pair(*i,eventRecord().incoming().second),clist);
  }

  if ( eventRecord().incoming().first->coloured() && eventRecord().incoming().second->coloured() ) {
    addCandidates(eventRecord().incoming(),clist);
  }

}

void DipoleShowerHandler::performSplitting(DipoleSplittingInfo& split) const {

  Ptr<DipoleSplittingKinematics>::tptr kinematics = split.splittingKinematics();
  kinematics->generateKinematics(split.emitter()->momentum(),
				 split.spectator()->momentum(),
				 split);

  split.splitEmitter(split.emitterData()->produceParticle(kinematics->lastEmitterMomentum()));
  split.splitSpectator(split.spectatorData()->produceParticle(kinematics->lastSpectatorMomentum()));
  split.emission(split.emissionData()->produceParticle(kinematics->lastEmissionMomentum()));

  // Setting resolution scales for the particles
  split.emission()->scale(sqr(split.lastPt()));
  split.splitEmitter()->scale(sqr(split.lastPt()));
  split.splitSpectator()->scale(split.spectator()->scale());

  PVector neighbours;
  if ( DipolePartonSplitter::colourConnected(split.emitter(),
					     eventRecord().incoming().first) &&
       split.emitter() != eventRecord().incoming().first )
    neighbours.push_back(eventRecord().incoming().first);
  if ( DipolePartonSplitter::colourConnected(split.emitter(),
					     eventRecord().incoming().second) &&
       split.emitter() != eventRecord().incoming().second )
    neighbours.push_back(eventRecord().incoming().second);
  for ( PList::const_iterator p = eventRecord().outgoing().begin(); 
	p != eventRecord().outgoing().end(); ++p ) {
    if ( *p == split.emitter() )
      continue;
    if ( DipolePartonSplitter::colourConnected(split.emitter(),*p) )
      neighbours.push_back(*p);
  }
  assert(neighbours.size() == 1 || neighbours.size() == 2 );
  if ( neighbours.size() == 2 ) {
    if ( UseRandom::rnd() < 0.5 )
      swap(neighbours[0],neighbours[1]);
  }

  DipolePartonSplitter::split(split.emitter(),split.splitEmitter(),split.emission(),
			      neighbours.front(),split.index().initialStateEmitter(),false);
  DipolePartonSplitter::change(split.spectator(),split.splitSpectator(),
			       split.index().initialStateSpectator(),false);
}

Energy DipoleShowerHandler::nextSubleadingSplitting(Energy hardPt,
						    Energy optHardPt, Energy optCutoff,
						    const bool decay) {

  list<DipoleSplittingInfo> candidates;
  getCandidates(candidates);

  hardScalesSubleading(candidates,hardPt);
  for ( list<DipoleSplittingInfo>::iterator cand = candidates.begin();
   	cand != candidates.end(); cand++ ) {
    cand->scale(hardPt);
  }


  list<DipoleSplittingInfo>::iterator split = candidates.end();

  // Winner of all dipoles
  DipoleSplittingInfo winner;
  // Winner for the current iteration of the for loop
  DipoleSplittingInfo candWinner;
  Energy winnerScale = 0.0*GeV;

  Energy nextScale = 0.0*GeV;

  for ( list<DipoleSplittingInfo>::iterator cand = candidates.begin();
	cand != candidates.end(); cand++ ) {
    nextScale = getWinner(candWinner,
			  cand->index(),
			  cand->emitterX(),cand->spectatorX(),
			  make_pair(true,false),
			  cand->emitter(),cand->spectator(),
			  hardPt,
			  optHardPt,
			  optCutoff);
    if ( nextScale > winnerScale ) {
      winnerScale = nextScale;
      winner = candWinner;
      split = cand;
      winnerIndex = winningKernelIndex;//check
    }
  }

  if ( split == candidates.end() )
    return ZERO;

  if ( decay )
    winner.isDecayProc( true );
  
  split->fill(winner);

  performSplitting(*split);
  eventRecord().update(*split);

  for ( list<DipoleSplittingInfo>::iterator dip = candidates.begin();
	dip != candidates.end(); ++dip ) {
    if ( dip == split )
      continue;
    dip->emission(split->emission());
    if ( dip->emitter() == split->emitter() ) {
      dip->splitEmitter(split->splitEmitter());
    } else {
      dip->splitEmitter(dip->emitter());
    }
    if ( dip->spectator() == split->spectator() ) {
      dip->splitSpectator(split->splitSpectator());
    } else {
      dip->splitSpectator(dip->spectator());
    }
  }

  // Update the ShowerHandler of the splitting reweight.
  if ( doSubleadingNc ) {
    theSplittingReweight->updateCurrentHandler();
  }

  return split->lastPt();

}


Energy DipoleShowerHandler::getWinner(DipoleSplittingInfo& winner,
                                      const Dipole& dip,
                                      pair<bool,bool> conf,
                                      Energy optHardPt,
                                      Energy optCutoff) {
  return
  getWinner(winner,dip.index(conf),
            dip.emitterX(conf),dip.spectatorX(conf),
            conf,dip.emitter(conf),dip.spectator(conf),
            dip.emitterScale(conf),optHardPt,optCutoff);
}

Energy DipoleShowerHandler::getWinner(SubleadingSplittingInfo& winner,
                                      Energy optHardPt,
                                      Energy optCutoff) {
  return
  getWinner(winner,winner.index(),
            winner.emitterX(),winner.spectatorX(),
            winner.configuration(),
            winner.emitter(),winner.spectator(),
            winner.startScale(),optHardPt,optCutoff);
}

Energy DipoleShowerHandler::getWinner(DipoleSplittingInfo& winner,
                                      const DipoleIndex& index,
                                      double emitterX, double spectatorX,
                                      pair<bool,bool> conf,
                                      tPPtr emitter, tPPtr spectator,
                                      Energy startScale,
                                      Energy optHardPt,
                                      Energy optCutoff) {
  
  if ( !index.initialStateEmitter() &&
      !doFSR() ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }
  
  if ( index.initialStateEmitter() &&
      !doISR() ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }

  if ( index.incomingDecaySpectator()
       && !doFSR() ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }
  
  // Currently do not split IF dipoles so
  // don't evaluate them in order to avoid
  // exceptions in the log
  if ( index.incomingDecayEmitter() ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }
  
  DipoleSplittingInfo candidate;
  candidate.index(index);
  candidate.configuration(conf);
  candidate.emitterX(emitterX);
  candidate.spectatorX(spectatorX);
  candidate.emitter(emitter);
  candidate.spectator(spectator);

  if ( generators().find(candidate.index()) == generators().end() )
  getGenerators(candidate.index(),theSplittingReweight);
  
  //
  // NOTE -- needs proper fixing at some point
  //
  // For some very strange reason, equal_range gives back
  // key ranges it hasn't been asked for. This particularly
  // happens e.g. for FI dipoles of the same kind, but different
  // PDF (hard vs MPI PDF). I can't see a reason for this,
  // as DipoleIndex properly implements comparison for equality
  // and (lexicographic) ordering; for the time being, we
  // use equal_range, extented by an explicit check for wether
  // the key is indeed what we wanted. See line after (*) comment
  // below.
  //
  // SW - Update 04/01/2016: Note - This caused a bug for me as I did not
  // include equality checks on the decay booleans in the == definition
  
  pair<GeneratorMap::iterator,GeneratorMap::iterator> gens
  = generators().equal_range(candidate.index());
  
  Energy winnerScale = 0.0*GeV;
  GeneratorMap::iterator winnerGen = generators().end();
  
  for ( GeneratorMap::iterator gen = gens.first; gen != gens.second; ++gen ) {

    if ( doPartialUnweighting )
      gen->second->doPartialUnweighting(referenceWeight);

    // (*) see NOTE above
    if ( !(gen->first == candidate.index()) )
    continue;
    
    if ( startScale <= gen->second->splittingKinematics()->IRCutoff() )
    continue;
    
    Energy dScale =
    gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),
                                                    spectator->momentum());

    // in very exceptional cases happening in DIS
    if ( std::isnan( double(dScale/MeV) ) )
    throw RedoShower();
    
    candidate.scale(dScale);
    
      // Calculate the mass of the recoil system
      // for decay dipoles
    if ( candidate.index().incomingDecaySpectator() || candidate.index().incomingDecayEmitter() ) {
      Energy recoilMass = gen->second->splittingKinematics()->recoilMassKin(emitter->momentum(),
                                                                            spectator->momentum());
      candidate.recoilMass(recoilMass);
    }

    // Store emitter and spectator masses, needed in kinematics
    if ( candidate.index().emitterData()->mass() != ZERO ) {
      if ( !candidate.index().offShellEmitter() )
	candidate.emitterMass( emitter->nominalMass() );  
      else
	candidate.emitterMass( emitter->mass() );  
    }
    
    if ( candidate.index().spectatorData()->mass() != ZERO ) {
      if ( !candidate.index().offShellSpectator() )
	candidate.spectatorMass( spectator->nominalMass() );  
      else
	candidate.spectatorMass( spectator->mass() );  
    }
    
    
    candidate.continuesEvolving();
    Energy hardScale = evolutionOrdering()->maxPt(startScale,candidate,*(gen->second->splittingKernel()));
    
    Energy maxPossible =
    gen->second->splittingKinematics()->ptMax(candidate.scale(),
                                              candidate.emitterX(), candidate.spectatorX(),
                                              candidate,
                                              *gen->second->splittingKernel());
    
    Energy ircutoff =
    optCutoff < gen->second->splittingKinematics()->IRCutoff() ?
    gen->second->splittingKinematics()->IRCutoff() :
    optCutoff;
    
    if ( maxPossible <= ircutoff ) {
      continue;
    }
    if ( maxPossible >= hardScale ){
      candidate.hardPt(hardScale);
    }
    else {
      hardScale = maxPossible;
      candidate.hardPt(maxPossible);
    }

    gen->second->generate(candidate,currentWeights(),optHardPt,optCutoff);
    Energy nextScale = evolutionOrdering()->evolutionScale(
                      gen->second->lastSplitting(),*(gen->second->splittingKernel()));
    
    if ( nextScale > winnerScale ) {
      winner.fill(candidate);
      gen->second->completeSplitting(winner);
      winnerGen = gen;
      winnerScale = nextScale;
      if ( continueSubleadingNc() )
	winningKernelIndex = kernelIndex+1;//check
    }
    
    if ( continueSubleadingNc() ) {
      kernelIndex++;//check
      scales.push_back(nextScale);//check
      theWeightsVector.push_back(gen->second->splittingWeightVector());
    }

    reweight(reweight() * gen->second->splittingWeight());
    
  }
  
  if ( winnerGen == generators().end() ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }
  
  if ( winner.stoppedEvolving() )
  return 0.0*GeV;
  
  return winnerScale;
  
}

void DipoleShowerHandler::doCascade(unsigned int& emDone,
                                    Energy optHardPt,
                                    Energy optCutoff,
                                    const bool decay) {
  
  if ( nEmissions )
    if ( emDone == nEmissions )
      return;

  if ( doSubleadingNc ) {
    unsigned int subEmDone = 0;
    // Set the starting scale
    Energy hardPt = muPt;

    double wref = referenceWeight;
    while ( subEmDone < subleadingNcEmissionsLimit && hardPt != ZERO && continueSubleadingNc() ) {
      // Clear out the weights from the earlier step
      theWeightsVector.clear();
      kernelIndex = 0;//check
      scales.clear();//check

      hardPt = nextSubleadingSplitting( hardPt, optHardPt, optCutoff, decay );

      // Partial unweighting
      if ( doPartialUnweightingAtEmission ) {
	const double w = reweight();
	if ( abs(w) < wref ) {
	  if ( abs(w)/wref < UseRandom::rnd() ) {
	    // Set weight to zero and end this event
	    reweight(0.0);
	    return;
	  } else
	    reweight( wref*w/abs(w) );
	}
	// Update the reference weight after emission
	wref *= referenceWeight;
      }

      // When the winning scale is larger than the cutoff
      // remove the added weights that are under the winning scale
      if ( hardPt != ZERO ) {
        Energy maxq = 0.0*GeV;
#ifndef NDEBUG
	size_t iwinner = theWeightsVector.size();//check
#endif
	for ( size_t i = 0; i < theWeightsVector.size(); i++ ) {
	  if ( theWeightsVector[i].size() > 0 ) {
	    // get<2> is true for an accept step.
	    if ( std::get<2>(theWeightsVector[i].back()) 
		 && std::get<0>(theWeightsVector[i].back()) > maxq) {
	      maxq = std::get<0>(theWeightsVector[i].back());
#ifndef NDEBUG
	      iwinner = i;//check
#endif
	    }
	  }
	}

        assert(winnerIndex-1 == iwinner);//check

	double correctionWeight = 1.0;
	for ( size_t i = 0; i < theWeightsVector.size(); i++ ) {
	  for ( size_t j = 0; j < theWeightsVector[i].size(); j++ ) {
	    if ( std::get<0>(theWeightsVector[i][j]) < maxq )
	      correctionWeight *= std::get<1>(theWeightsVector[i][j]);
	  }
	}
	reweight(reweight()/correctionWeight);
      }

      // Increment the number of subleading Nc emissions done
      subEmDone++;
      // Stop if the limit of emissions is reached
      if ( nEmissions )
	if ( ++emDone == nEmissions )
	  return;
    }

    // Subleading shower done, prepare chains for the standard
    // dipole shower
    eventRecord().prepareChainsSubleading( decay );
    // Set scales
    for ( list<DipoleChain>::iterator ch = eventRecord().chains().begin();
	  ch != eventRecord().chains().end(); ch++ ) {
      for ( list<Dipole>::iterator dp = ch->dipoles().begin();
	    dp != ch->dipoles().end(); dp++ ) {
	dp->emitterScale(make_pair(true,false),hardPt);
	dp->emitterScale(make_pair(false,true),hardPt);
      }
    }
    
  }


  DipoleSplittingInfo winner;
  DipoleSplittingInfo dipoleWinner;
  
  
  while ( eventRecord().haveChain() ) {

    if ( verbosity > 2 ) {
      generator()->log() << "DipoleShowerHandler selecting splittings for the chain:\n"
      << eventRecord().currentChain() << flush;
    }
    
    list<Dipole>::iterator winnerDip = eventRecord().currentChain().dipoles().end();
    Energy winnerScale = 0.0*GeV;
    Energy nextLeftScale = 0.0*GeV;
    Energy nextRightScale = 0.0*GeV;
    for ( list<Dipole>::iterator dip = eventRecord().currentChain().dipoles().begin();
         dip != eventRecord().currentChain().dipoles().end(); ++dip ) {
      
      nextLeftScale = getWinner(dipoleWinner,*dip,{true,false},optHardPt,optCutoff);
      if ( nextLeftScale > winnerScale ) {
        winnerScale = nextLeftScale;
        winner = dipoleWinner;
        winnerDip = dip;
      }
      
      nextRightScale = getWinner(dipoleWinner,*dip,{false,true},optHardPt,optCutoff);
      if ( nextRightScale > winnerScale ) {
        winnerScale = nextRightScale;
        winner = dipoleWinner;
        winnerDip = dip;
      }
      
      if ( evolutionOrdering()->independentDipoles() ) {
        Energy dipScale = max(nextLeftScale,nextRightScale);
        if ( dip->leftScale() > dipScale )
        dip->leftScale(dipScale);
        if ( dip->rightScale() > dipScale )
        dip->rightScale(dipScale);
      }
    }
    
    if ( verbosity > 1 ) {
      if ( winnerDip != eventRecord().currentChain().dipoles().end() )
      generator()->log() << "DipoleShowerHandler selected the splitting:\n"
			   << winner << " for the dipole\n"
			   << (*winnerDip) << flush;
      else
      generator()->log() << "DipoleShowerHandler could not select a splitting above the IR cutoff\n"
			   << flush;
    }
    
      // pop the chain if no dipole did radiate
    if ( winnerDip == eventRecord().currentChain().dipoles().end() ) {
      eventRecord().popChain();
      if ( theEventReweight && eventRecord().chains().empty() )
      if ( (theEventReweight->firstInteraction() && firstInteraction()) ||
          (theEventReweight->secondaryInteractions() && !firstInteraction()) ) {
        double w = theEventReweight->weightCascade(eventRecord().incoming(),
                                                   eventRecord().outgoing(),
                                                   eventRecord().hard(),theGlobalAlphaS);
        reweight(reweight()*w);
      }
      continue;
    }
    
      // otherwise perform the splitting
      // but first see if the emission would produce a configuration in the ME region.
    if (   theMergingHelper
	&& eventHandler()->currentCollision()
	&& !decay
	&& firstInteraction() ) {
      if (theMergingHelper->maxLegs()>eventRecord().outgoing().size()+
                                      eventRecord().hard().size()
                                      +2){//incoming
        
	if (theMergingHelper->mergingScale()<winnerScale && 
            theMergingHelper->emissionProbability() < UseRandom::rnd()) {
            
          theMergingHelper->setEmissionProbability(0.);

          const bool transparent=true;
          if (transparent) {
            pair<list<Dipole>::iterator,list<Dipole>::iterator> tmpchildren;
            DipoleSplittingInfo tmpwinner=winner;
            DipoleChain* tmpfirstChain = nullptr;
            DipoleChain* tmpsecondChain = nullptr;
            
            auto New=eventRecord().tmpsplit(winnerDip,tmpwinner,
                                            tmpchildren,tmpfirstChain,
                                            tmpsecondChain);
            
            
            if (theMergingHelper->matrixElementRegion(New.first,
                                                      New.second,
                                                      winnerScale,
                                                      theMergingHelper->mergingScale())) {
              optHardPt=winnerScale;
              continue;
            }
          }else{
            optHardPt=winnerScale;
            continue;
            
          }
        }
      }
    }
    if(theMergingHelper&&firstInteraction())
       optHardPt=ZERO;  
  
   
    
    didRadiate = true;
    
    eventRecord().isMCatNLOSEvent(false);
    eventRecord().isMCatNLOHEvent(false);
    
    pair<list<Dipole>::iterator,list<Dipole>::iterator> children;
    
    DipoleChain* firstChain = nullptr;
    DipoleChain* secondChain = nullptr;
    
    // Generate the azimuthal angle 
    if ( spinCorrelations() ) 
      vertexRecord().generatePhi(winner,*winnerDip);

    if ( decay )
      winner.isDecayProc( true );
    
      // Note: the dipoles are updated in eventRecord().split(....) after the splitting,
      // hence the entire cascade is handled in doCascade
      // The dipole scales are updated in dip->split(....)
    
    if ( decay )
    winner.isDecayProc( true );

    eventRecord().split(winnerDip,winner,children,firstChain,secondChain);

    // Update the vertex record following the splitting
    if ( spinCorrelations() )
      vertexRecord().update(winner);

    assert(firstChain && secondChain);
    evolutionOrdering()->setEvolutionScale(winnerScale,winner,*firstChain,children);
    if ( !secondChain->dipoles().empty() )
    evolutionOrdering()->setEvolutionScale(winnerScale,winner,*secondChain,children);
    
    if ( verbosity > 1 ) {
      generator()->log() << "DipoleShowerHandler did split the last selected dipole into:\n"
      << (*children.first) << (*children.second) << flush;
    }
    
    if ( verbosity > 2 ) {
      generator()->log() << "After splitting the last selected dipole, "
      << "DipoleShowerHandler encountered the following chains:\n"
      << (*firstChain) << (*secondChain) << flush;
    }
    
    if ( theEventReweight )
    if ( (theEventReweight->firstInteraction() && firstInteraction()) ||
        (theEventReweight->secondaryInteractions() && !firstInteraction()) ) {
      double w = theEventReweight->weight(eventRecord().incoming(),
                                          eventRecord().outgoing(),
                                          eventRecord().hard(),theGlobalAlphaS);
      reweight(reweight()*w);
    }
    
    if ( nEmissions )
    if ( ++emDone == nEmissions )
    return;
  }
  
}

bool DipoleShowerHandler::realign() {
  
  if ( !didRadiate && !intrinsicPtGenerator )
  return false;
  
  if ( eventRecord().incoming().first->coloured() ||
      eventRecord().incoming().second->coloured() ) {
    
    if ( eventRecord().incoming().first->momentum().perp2()/GeV2 < 1e-10 &&
        eventRecord().incoming().second->momentum().perp2()/GeV2 < 1e-10 )
    return false;
    
    pair<Lorentz5Momentum,Lorentz5Momentum> inMomenta
    (eventRecord().incoming().first->momentum(),
     eventRecord().incoming().second->momentum());
    
    LorentzRotation transform((inMomenta.first+inMomenta.second).findBoostToCM());
    
    Axis dir = (transform * inMomenta.first).vect().unit();
    Axis rot (-dir.y(),dir.x(),0);
    double theta = dir.theta();
    
    if ( lastParticles().first->momentum().z() < ZERO )
    theta = -theta;
    
    transform.rotate(-theta,rot);
    
    inMomenta.first = transform*inMomenta.first;
    inMomenta.second = transform*inMomenta.second;
    
    assert(inMomenta.first.z() > ZERO &&
           inMomenta.second.z() < ZERO);
    
    Energy2 sHat =
    (eventRecord().incoming().first->momentum() +
     eventRecord().incoming().second->momentum()).m2();
    
    pair<Energy,Energy> masses(eventRecord().incoming().first->mass(),
                               eventRecord().incoming().second->mass());
    pair<Energy,Energy> qs;
    
    if ( !eventRecord().incoming().first->coloured() ) {
      assert(masses.second == ZERO);
      qs.first = eventRecord().incoming().first->momentum().z();
      qs.second = (sHat-sqr(masses.first))/(2.*(qs.first+sqrt(sqr(masses.first)+sqr(qs.first))));
    } else if ( !eventRecord().incoming().second->coloured() ) {
      assert(masses.first == ZERO);
      qs.second = eventRecord().incoming().second->momentum().z();
      qs.first = (sHat-sqr(masses.second))/(2.*(qs.second+sqrt(sqr(masses.second)+sqr(qs.second))));
    } else {
      assert(masses.first == ZERO && masses.second == ZERO);
      if ( realignmentScheme == 0 ) {
        double yX = eventRecord().pX().rapidity();
        double yInt = (transform*eventRecord().pX()).rapidity();
        double dy = yX-yInt;
        qs.first = (sqrt(sHat)/2.)*exp(dy);
        qs.second = (sqrt(sHat)/2.)*exp(-dy);
      } else if ( realignmentScheme == 1 ) {
        Energy sS = sqrt((lastParticles().first->momentum() +
                          lastParticles().second->momentum()).m2());
        qs.first = eventRecord().fractions().first * sS / 2.;
        qs.second = eventRecord().fractions().second * sS / 2.;
      }
    }
    
    double beta =
    (qs.first-qs.second) /
    ( sqrt(sqr(masses.first)+sqr(qs.first)) +
     sqrt(sqr(masses.second)+sqr(qs.second)) );
    transform.boostZ(beta);
    
    Lorentz5Momentum tmp;
    
    if ( eventRecord().incoming().first->coloured() ) {
      tmp = eventRecord().incoming().first->momentum();
      tmp = transform * tmp;
      eventRecord().incoming().first->set5Momentum(tmp);
    }
    if ( eventRecord().incoming().second->coloured() ) {
      tmp = eventRecord().incoming().second->momentum();
      tmp = transform * tmp;
      eventRecord().incoming().second->set5Momentum(tmp);
    }
    eventRecord().transform(transform);
    return true;
    
  }
  
  return false;
  
}

void DipoleShowerHandler::resetAlphaS(Ptr<AlphaSBase>::tptr as) {
  
  for ( auto & k : kernels) {
    if ( !k->alphaS() )
    k->alphaS(as);
    k->renormalizationScaleFreeze(theRenormalizationScaleFreeze);
    k->factorizationScaleFreeze(theFactorizationScaleFreeze);
  }
  
    // clear the generators to be rebuild
    // actually, there shouldn't be any generators
    // when this happens.
  generators().clear();
  
}

void DipoleShowerHandler::resetReweight(Ptr<DipoleSplittingReweight>::tptr rw) {
  for ( auto &  g : generators() )
  g.second->splittingReweight(rw);
}

void DipoleShowerHandler::getGenerators(const DipoleIndex& ind,
                                        Ptr<DipoleSplittingReweight>::tptr rw) {
  
  bool gotone = false;
  
  for ( auto &  k : kernels ) {
    if ( k->canHandle(ind) ) {
      
      if ( verbosity > 0 ) {
        generator()->log() << "DipoleShowerHandler encountered the dipole configuration\n"
        << ind << " in event number "
        << eventHandler()->currentEvent()->number()
        << "\nwhich can be handled by the splitting kernel '"
        << k->name() << "'.\n" << flush;
      }
      
      gotone = true;
      Ptr<DipoleSplittingGenerator>::ptr nGenerator =
      new_ptr(DipoleSplittingGenerator());
      nGenerator->doCompensate(theDoCompensate);
      nGenerator->splittingKernel(k);
      if ( renormalizationScaleFactor() != 1. )
      nGenerator->splittingKernel()->renormalizationScaleFactor(renormalizationScaleFactor());
      if ( factorizationScaleFactor() != 1. )
      nGenerator->splittingKernel()->factorizationScaleFactor(factorizationScaleFactor());
      if ( !nGenerator->splittingReweight() )
      nGenerator->splittingReweight(rw);
      nGenerator->splittingKernel()->freezeGrid(theFreezeGrid);
      nGenerator->splittingKernel()->detuning(theDetuning);
      
      GeneratorMap::const_iterator equivalent = generators().end();
      
      for ( GeneratorMap::const_iterator eq = generators().begin();
           eq != generators().end(); ++eq ) {
        if ( !eq->second->wrapping() )
        if ( k->canHandleEquivalent(ind,*(eq->second->splittingKernel()),eq->first) ) {
          
          equivalent = eq;
          
          if ( verbosity > 0 ) {
            generator()->log() << "The dipole configuration "
            << ind
            << " can equivalently be handled by the existing\n"
            << "generator for configuration "
            << eq->first << " using the kernel '"
            << eq->second->splittingKernel()->name()
            << "'\n" << flush;
          }
          
          break;
          
        }
      }
      
      if ( equivalent != generators().end() ) {
        nGenerator->wrap(equivalent->second);
      }
      
      DipoleSplittingInfo dummy;
      dummy.index(ind);
      nGenerator->prepare(dummy);
      
      generators().insert({ind,nGenerator});
      
    }
  }
  
  if ( !gotone ) {
    throw Exception()
      << "DipoleShowerHandler could not "
      << "find a splitting kernel which is able "
      << "to handle splittings off the dipole "
      << ind << ".\n"
      << "Please check the input files."
      << Exception::runerror;
  }
  
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void DipoleShowerHandler::doinit() {
  ShowerHandler::doinit();
  if ( theGlobalAlphaS )
  resetAlphaS(theGlobalAlphaS);
  // copy off-shell particle ids before showering from input vector to the 
  // set used in the simulation
  if ( theColouredOffShellInShower.empty() ) {
    for(unsigned int ix=0;ix<theInputColouredOffShellInShower.size();++ix)
      theColouredOffShellInShower.insert(abs(theInputColouredOffShellInShower[ix]));
  }
  // work out which shower phase space to use for the matching
  bool zChoice0 = false;
  bool zChoice1 = false;
  size_t zChoiceOther = false;
  for ( auto & k : kernels) {
    if ( k->splittingKinematics()->openZBoundaries() == 0 )
      zChoice0 = true;
    else if ( k->splittingKinematics()->openZBoundaries() == 1 )
      zChoice1 = true;
    else
      zChoiceOther = true;
    // either inconsistent or other option which cannot be handled by the matching
    if ( zChoice0 && zChoice1 ) {
      zChoiceOther = true; break;
    }
  }
  if ( zChoiceOther )
    theZBoundaries = 2;
  else if ( zChoice1 )
    theZBoundaries = 1;
  else if ( zChoice0 )
    theZBoundaries = 0;
}

void DipoleShowerHandler::dofinish() {
  ShowerHandler::dofinish();
}

void DipoleShowerHandler::doinitrun() {
  ShowerHandler::doinitrun();
}

void DipoleShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << kernels << theEvolutionOrdering
     << constituentReshuffler << intrinsicPtGenerator
     << theGlobalAlphaS << chainOrderVetoScales
     << nEmissions << discardNoEmissions << firstMCatNLOEmission
     << thePowhegDecayEmission
    //<< theAnalyseSpinCorrelations
     << doSubleadingNc << subleadingNcEmissionsLimit
     << densityOperatorEvolution << ounit(densityOperatorCutoff,GeV2)
     << doPartialUnweightingAtEmission
     << doPartialUnweighting << referenceWeight
     << cmecReweightFactor << negCMECScaling
     << realignmentScheme << verbosity << printEvent
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV)
     << theShowerApproximation
     << theDoCompensate << theFreezeGrid << theDetuning
     << theEventReweight << theSplittingReweight << ounit(maxPt,GeV)
     << ounit(muPt,GeV)<< theMergingHelper << theColouredOffShellInShower
     << theInputColouredOffShellInShower
      << theZBoundaries;
}

void DipoleShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> kernels >> theEvolutionOrdering
     >> constituentReshuffler >> intrinsicPtGenerator
     >> theGlobalAlphaS >> chainOrderVetoScales
     >> nEmissions >> discardNoEmissions >> firstMCatNLOEmission
     >> thePowhegDecayEmission
    //>> theAnalyseSpinCorrelations
     >> doSubleadingNc >> subleadingNcEmissionsLimit
     >> densityOperatorEvolution >> iunit(densityOperatorCutoff,GeV2)
     >> doPartialUnweightingAtEmission
     >> doPartialUnweighting >> referenceWeight
     >> cmecReweightFactor >> negCMECScaling
     >> realignmentScheme >> verbosity >> printEvent
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV)
     >> theShowerApproximation
     >> theDoCompensate >> theFreezeGrid >> theDetuning
     >> theEventReweight >> theSplittingReweight >> iunit(maxPt,GeV)
     >> iunit(muPt,GeV)>>theMergingHelper >> theColouredOffShellInShower
     >> theInputColouredOffShellInShower
      >> theZBoundaries;
}

ClassDescription<DipoleShowerHandler> DipoleShowerHandler::initDipoleShowerHandler;
  // Definition of the static class description member.

void DipoleShowerHandler::Init() {
  
  static ClassDocumentation<DipoleShowerHandler> documentation
  ("The DipoleShowerHandler class manages the showering using "
   "the dipole shower algorithm.",
   "The shower evolution was performed using the algorithm described in "
   "\\cite{Platzer:2009jq} and \\cite{Platzer:2011bc}.",
   "%\\cite{Platzer:2009jq}\n"
   "\\bibitem{Platzer:2009jq}\n"
   "S.~Platzer and S.~Gieseke,\n"
   "``Coherent Parton Showers with Local Recoils,''\n"
   "  JHEP {\\bf 1101}, 024 (2011)\n"
   "arXiv:0909.5593 [hep-ph].\n"
   "%%CITATION = ARXIV:0909.5593;%%\n"
   "%\\cite{Platzer:2011bc}\n"
   "\\bibitem{Platzer:2011bc}\n"
   "S.~Platzer and S.~Gieseke,\n"
   "``Dipole Showers and Automated NLO Matching in Herwig,''\n"
   "arXiv:1109.6256 [hep-ph].\n"
   "%%CITATION = ARXIV:1109.6256;%%");
  
  static RefVector<DipoleShowerHandler,DipoleSplittingKernel> interfaceKernels
  ("Kernels",
   "Set the splitting kernels to be used by the dipole shower.",
   &DipoleShowerHandler::kernels, -1, false, false, true, false, false);
  
  
  static Reference<DipoleShowerHandler,DipoleEvolutionOrdering> interfaceEvolutionOrdering
  ("EvolutionOrdering",
   "Set the evolution ordering to be used.",
   &DipoleShowerHandler::theEvolutionOrdering, false, false, true, false, false);
  
  
  static Reference<DipoleShowerHandler,ConstituentReshuffler> interfaceConstituentReshuffler
  ("ConstituentReshuffler",
   "The object to be used to reshuffle partons to their constitutent mass shells.",
   &DipoleShowerHandler::constituentReshuffler, false, false, true, true, false);
  
  
  static Reference<DipoleShowerHandler,IntrinsicPtGenerator> interfaceIntrinsicPtGenerator
  ("IntrinsicPtGenerator",
   "Set the object in charge to generate intrinsic pt for incoming partons.",
   &DipoleShowerHandler::intrinsicPtGenerator, false, false, true, true, false);
  
  static Reference<DipoleShowerHandler,AlphaSBase> interfaceGlobalAlphaS
    ("GlobalAlphaS",
     "Set a global strong coupling for all splitting kernels.",
     &DipoleShowerHandler::theGlobalAlphaS, false, false, true, true, false);
// Start: Trying to add interface for the subleading Nc
  static Switch<DipoleShowerHandler,bool> interfaceDoSubleadingNc
    ("DoSubleadingNc",
     "Switch on or off subleading Nc corrections.",
     &DipoleShowerHandler::doSubleadingNc, true, false, false);
  static SwitchOption interfaceDoSubleadingNcOn
    (interfaceDoSubleadingNc,
     "On",
     "Switch on subleading Nc corrections.",
     true);
  static SwitchOption interfaceDoSubleadingNcOff
    (interfaceDoSubleadingNc,
     "Off",
     "Switch off subleading Nc corrections.",
     false);
  // Limit for how many subleading Nc emissions should be calculated
  static Parameter<DipoleShowerHandler,size_t> interfaceSubleadingNcEmissionsLimit
    ("SubleadingNcEmissionsLimit",
     "Number of emissions to calculate subleading Nc corrections for.",
     &DipoleShowerHandler::subleadingNcEmissionsLimit,0,0,0,
     false, false, Interface::lowerlim);
  static Parameter<DipoleShowerHandler,int> interfaceDensityOperatorEvolution
    ("DensityOperatorEvolution",
     "Scheme for evolving the density operator.",
     &DipoleShowerHandler::densityOperatorEvolution,0,0,0,
     false, false, Interface::lowerlim);
  static Parameter<DipoleShowerHandler,Energy2> interfaceDensityOperatorCutoff
    ("DensityOperatorCutoff",
     "Cutoff for momentum invariants for the density operator evolution.",
     &DipoleShowerHandler::densityOperatorCutoff,GeV2,1.0*GeV2,0.0*GeV2,0*GeV2,
     false, false, Interface::lowerlim);
  static Switch<DipoleShowerHandler,bool> interfaceDoPartialUnweightingAtEmission
    ("DoPartialUnweightingAtEmission",
     "Switch on or off partial unweighting at the emission level.",
     &DipoleShowerHandler::doPartialUnweightingAtEmission,true,false,false);
  static SwitchOption interfaceDoPartialUnweightingAtEmissionOn
    (interfaceDoPartialUnweightingAtEmission,
     "On",
     "Switch on partial unweighting.",
     true);
  static SwitchOption interfaceDoPartialUnweightingAtEmissionOff
    (interfaceDoPartialUnweightingAtEmission,
     "Off",
     "Switch off partial unweighting.",
     false);
  static Switch<DipoleShowerHandler,bool> interfaceDoPartialUnweighting
    ("DoPartialUnweighting",
     "Switch on or off partial unweighting at the dipole splitting level.",
     &DipoleShowerHandler::doPartialUnweighting,true,false,false);
  static SwitchOption interfaceDoPartialUnweightingOn
    (interfaceDoPartialUnweighting,
     "On",
     "Switch on partial unweighting.",
     true);
  static SwitchOption interfaceDoPartialUnweightingOff
    (interfaceDoPartialUnweighting,
     "Off",
     "Switch off partial unweighting.",
     false);
  static Parameter<DipoleShowerHandler,double> interfaceReferenceWeight
    ("ReferenceWeight",
     "Reference weight for the partial unweighting.",
     &DipoleShowerHandler::referenceWeight,0.1,0.0,0,
     false, false, Interface::lowerlim);
  static Parameter<DipoleShowerHandler,double> interfaceCMECReweightFactor
    ("CMECReweightFactor",
     "Factor used in the reweighting algorithm.",
     &DipoleShowerHandler::cmecReweightFactor,1.0,0.0,0,
     false, false, Interface::lowerlim);
  static Parameter<DipoleShowerHandler,double> interfaceNegCMECScaling
    ("NegCMECScaling",
     "Scaling factor for the negative colour matrix element corrections (CMECs).",
     &DipoleShowerHandler::negCMECScaling,0.0,0.0,0,
     false, false, Interface::lowerlim);
  static Switch<DipoleShowerHandler,int> interfaceRealignmentScheme
  ("RealignmentScheme",
   "The realignment scheme to use.",
   &DipoleShowerHandler::realignmentScheme, 0, false, false);
  static SwitchOption interfaceRealignmentSchemePreserveRapidity
  (interfaceRealignmentScheme,
   "PreserveRapidity",
   "Preserve the rapidity of non-coloured outgoing system.",
   0);
  static SwitchOption interfaceRealignmentSchemeEvolutionFractions
  (interfaceRealignmentScheme,
   "EvolutionFractions",
   "Use momentum fractions as generated by the evolution.",
   1);
  static SwitchOption interfaceRealignmentSchemeCollisionFrame
  (interfaceRealignmentScheme,
   "CollisionFrame",
   "Determine realignment from collision frame.",
   2);
  
  
  static Switch<DipoleShowerHandler,bool> interfaceChainOrderVetoScales
  ("ChainOrderVetoScales",
   "[experimental] Switch the chain ordering for veto scales on or off.",
   &DipoleShowerHandler::chainOrderVetoScales, true, false, false);
  static SwitchOption interfaceChainOrderVetoScalesYes
  (interfaceChainOrderVetoScales,
   "Yes",
   "Switch on chain ordering for veto scales.",
   true);
  static SwitchOption interfaceChainOrderVetoScalesNo
  (interfaceChainOrderVetoScales,
   "No",
   "Switch off chain ordering for veto scales.",
   false);
  
  interfaceChainOrderVetoScales.rank(-1);
  
  
  static Parameter<DipoleShowerHandler,unsigned int> interfaceNEmissions
  ("NEmissions",
   "[debug option] Limit the number of emissions to be generated. Zero does not limit the number of emissions.",
   &DipoleShowerHandler::nEmissions, 0, 0, 0,
   false, false, Interface::lowerlim);
  
  interfaceNEmissions.rank(-1);
  
  
  static Switch<DipoleShowerHandler,bool> interfaceDiscardNoEmissions
  ("DiscardNoEmissions",
   "[debug option] Discard events without radiation.",
   &DipoleShowerHandler::discardNoEmissions, false, false, false);
  static SwitchOption interfaceDiscardNoEmissionsYes
  (interfaceDiscardNoEmissions,
   "Yes",
   "Discard events without radiation.",
   true);
  static SwitchOption interfaceDiscardNoEmissionsNo
  (interfaceDiscardNoEmissions,
   "No",
   "Do not discard events without radiation.",
   false);
  
  interfaceDiscardNoEmissions.rank(-1);
  
  static Switch<DipoleShowerHandler,bool> interfaceFirstMCatNLOEmission
  ("FirstMCatNLOEmission",
   "[debug option] Only perform the first MC@NLO emission.",
   &DipoleShowerHandler::firstMCatNLOEmission, false, false, false);
  static SwitchOption interfaceFirstMCatNLOEmissionYes
  (interfaceFirstMCatNLOEmission,
   "Yes",
   "Perform only the first MC@NLO emission.",
   true);
  static SwitchOption interfaceFirstMCatNLOEmissionNo
  (interfaceFirstMCatNLOEmission,
   "No",
   "Produce all emissions.",
   false);
  
  interfaceFirstMCatNLOEmission.rank(-1);
  
  
  static Parameter<DipoleShowerHandler,int> interfaceVerbosity
  ("Verbosity",
   "[debug option] Set the level of debug information provided.",
   &DipoleShowerHandler::verbosity, 0, 0, 0,
   false, false, Interface::lowerlim);
  
  interfaceVerbosity.rank(-1);
  
  
  static Parameter<DipoleShowerHandler,int> interfacePrintEvent
  ("PrintEvent",
   "[debug option] The number of events for which debugging information should be provided.",
   &DipoleShowerHandler::printEvent, 0, 0, 0,
   false, false, Interface::lowerlim);
  
  interfacePrintEvent.rank(-1);
  
  static Parameter<DipoleShowerHandler,Energy> interfaceRenormalizationScaleFreeze
  ("RenormalizationScaleFreeze",
   "The freezing scale for the renormalization scale.",
   &DipoleShowerHandler::theRenormalizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
   false, false, Interface::lowerlim);
  
  static Parameter<DipoleShowerHandler,Energy> interfaceFactorizationScaleFreeze
  ("FactorizationScaleFreeze",
   "The freezing scale for the factorization scale.",
   &DipoleShowerHandler::theFactorizationScaleFreeze, GeV, 2.0*GeV, 0.0*GeV, 0*GeV,
   false, false, Interface::lowerlim);
  
  static Switch<DipoleShowerHandler,bool> interfaceDoCompensate
  ("DoCompensate",
   "",
   &DipoleShowerHandler::theDoCompensate, false, false, false);
  static SwitchOption interfaceDoCompensateYes
  (interfaceDoCompensate,
   "Yes",
   "",
   true);
  static SwitchOption interfaceDoCompensateNo
  (interfaceDoCompensate,
   "No",
   "",
   false);
  
  static Parameter<DipoleShowerHandler,unsigned long> interfaceFreezeGrid
  ("FreezeGrid",
   "",
   &DipoleShowerHandler::theFreezeGrid, 500000, 1, 0,
   false, false, Interface::lowerlim);
  
  static Parameter<DipoleShowerHandler,double> interfaceDetuning
  ("Detuning",
   "A value to detune the overestimate kernel.",
   &DipoleShowerHandler::theDetuning, 1.0, 1.0, 0,
   false, false, Interface::lowerlim);
  
  static Reference<DipoleShowerHandler,DipoleEventReweight> interfaceEventReweight
  ("EventReweight",
   "",
   &DipoleShowerHandler::theEventReweight, false, false, true, true, false);
  
  static Reference<DipoleShowerHandler,DipoleSplittingReweight> interfaceSplittingReweight
  ("SplittingReweight",
   "Set the splitting reweight.",
   &DipoleShowerHandler::theSplittingReweight, false, false, true, true, false);
  
  
  static Switch<DipoleShowerHandler, bool> interfacePowhegDecayEmission
  ("PowhegDecayEmission",
   "Use Powheg style emission for the decays",
   &DipoleShowerHandler::thePowhegDecayEmission, true, false, false);
  
  static SwitchOption interfacePowhegDecayEmissionYes
  (interfacePowhegDecayEmission,"Yes","Powheg decay emission on", true);
  
  static SwitchOption interfacePowhegDecayEmissionNo
  (interfacePowhegDecayEmission,"No","Powheg decay emission off", false);

  static ParVector<DipoleShowerHandler,long> interfaceOffShellInShower
    ("OffShellInShower",
     "PDG codes of the coloured particles that can be off-shell in the process.",
     &DipoleShowerHandler::theInputColouredOffShellInShower, -1, 0l, -10000000l, 10000000l,
     false, false, Interface::limited);

  /*
    static Switch<DipoleShowerHandler, bool> interfaceAnalyseSpinCorrelations
    ("AnalyseSpinCorrelations",
    "Record the information required for the spin correlation analyis.",
    &DipoleShowerHandler::theAnalyseSpinCorrelations, false, false, false);
  
    static SwitchOption interfaceAnalyseSpinCorrelationsYes
    (interfaceAnalyseSpinCorrelations,"Yes","Record the information for analysing the spin correlations.", true);
  
    static SwitchOption interfaceAnalyseSpinCorrelationsNo
    (interfaceAnalyseSpinCorrelations,"No","Do not record extra information.", false);
  */

  

  
}

