// -*- C++ -*-
//
// DipoleShowerHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

// include theses to have complete types
#include "Herwig/Shower/Base/Evolver.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"

#include "Herwig/DipoleShower/Utility/DipolePartonSplitter.h"

#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

bool DipoleShowerHandler::firstWarn = true;

DipoleShowerHandler::DipoleShowerHandler() 
  : ShowerHandler(), chainOrderVetoScales(true),
    nEmissions(0), discardNoEmissions(false), firstMCatNLOEmission(false),
    doFSR(true), doISR(true), realignmentScheme(0),
    verbosity(0), printEvent(0), nTries(0), 
    didRadiate(false), didRealign(false),
    theRenormalizationScaleFreeze(1.*GeV), 
    theFactorizationScaleFreeze(2.*GeV),
    isMCatNLOSEvent(false),
    isMCatNLOHEvent(false), theDoCompensate(false),
    theFreezeGrid(500000), theDetuning(1.0),
    maxPt(ZERO), muPt(ZERO) {}

DipoleShowerHandler::~DipoleShowerHandler() {}

IBPtr DipoleShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleShowerHandler::fullclone() const {
  return new_ptr(*this);
}

tPPair DipoleShowerHandler::cascade(tSubProPtr sub, XCPtr,
				    Energy optHardPt, Energy optCutoff) {

  useMe();

  prepareCascade(sub);
  resetWeights();

  if ( !doFSR && ! doISR )
    return sub->incoming();

  eventRecord().clear();
  eventRecord().prepare(sub,dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr()),pdfs());

  if ( eventRecord().outgoing().empty() && !doISR )
    return sub->incoming();
  if ( !eventRecord().incoming().first->coloured() &&
       !eventRecord().incoming().second->coloured() &&
       !doFSR )
    return sub->incoming();

  nTries = 0;

  while ( true ) {

    try {

      didRadiate = false;
      didRealign = false;

      isMCatNLOSEvent = false;
      isMCatNLOHEvent = false;

      if ( eventRecord().xcombPtr() ) {

	Ptr<SubtractedME>::tptr subme =
	  dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(eventRecord().xcombPtr()->matrixElement());
	Ptr<MatchboxMEBase>::tptr me =
	  dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(eventRecord().xcombPtr()->matrixElement());
	Ptr<SubtractionDipole>::tptr dipme =  
	  dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(eventRecord().xcombPtr()->matrixElement());


	if ( subme ) {
	  if ( subme->showerApproximation() ) {
	    // don't do this for POWHEG-type corrections
	    if ( !subme->showerApproximation()->needsSplittingGenerator() ) {
	      theShowerApproximation = subme->showerApproximation();
	      if ( subme->realShowerSubtraction() )
		isMCatNLOHEvent = true;
	      else if ( subme->virtualShowerSubtraction() )
		isMCatNLOSEvent = true;
	    }
	  }
	} else if ( me ) {
	  if ( me->factory()->showerApproximation() ) {
	    if ( !me->factory()->showerApproximation()->needsSplittingGenerator() ) {
	      theShowerApproximation = me->factory()->showerApproximation();
	      isMCatNLOSEvent = true;
	    }
	  }
	}

	string error = "Inconsistent hard emission set-up in DipoleShowerHandler::cascade. ";      
	if (evolver()->hardEmissionMode()==1 || evolver()->hardEmissionMode()==3 )
	  throw Exception() << error
			    << "Cannot generate POWHEG corrections "
			    << "for particle decays using DipoleShowerHandler.  "
			    << "Check value of Evolver:HardEmissionMode."
			    << Exception::runerror; 
	if ( ( isMCatNLOSEvent || isMCatNLOHEvent ) &&  evolver()->hardEmissionMode()==2)
	  throw Exception() << error
			    << "Cannot generate POWHEG matching with MC@NLO shower approximation. "
			    << "Add 'set Evolver:HardEmissionMode 0' to input file."
			    << Exception::runerror;
	if (me && me->factory()->showerApproximation()){
	  if(me->factory()->showerApproximation()->needsTruncatedShower())
	    throw Exception() << error
			      << "No truncated shower needed with DipoleShowerHandler.  Add "
			      << "'set MEMatching:TruncatedShower No' to input file."
			      << Exception::runerror;
	  if (!( isMCatNLOSEvent || isMCatNLOHEvent ) &&
	      evolver()->hardEmissionMode()==0 && firstWarn){
	    firstWarn=false;
	    throw Exception() << error
			      << "Evolver:HardEmissionMode will be set to 'MatchboxPOWHEG'"
			      << Exception::warning;
	  }  
	}
	else if (subme && subme->factory()->showerApproximation()){
	  if(subme->factory()->showerApproximation()->needsTruncatedShower())
	    throw Exception() << error
			      << "No truncated shower needed with DipoleShowerHandler.  Add "
			      << "'set MEMatching:TruncatedShower No' to input file."
			      << Exception::runerror;
	  if (!( isMCatNLOSEvent || isMCatNLOHEvent ) &&
	      evolver()->hardEmissionMode()==0 && firstWarn){
	    firstWarn=false;
	    throw Exception() << error
			      << "Evolver:HardEmissionMode will be set to 'MatchboxPOWHEG'"
			      << Exception::warning;
	  }  
	}
	else if (dipme && evolver()->hardEmissionMode() == 0 && firstWarn){
	  firstWarn=false;
	  throw Exception() << error
			    << "Evolver:HardEmissionMode will be set to 'MatchboxPOWHEG'"
			    << Exception::warning;
	}
	else if (!dipme && evolver()->hardEmissionMode()==2 &&
		 ShowerHandler::currentHandler()->firstInteraction())
	  throw Exception() << error
			    << "POWHEG matching requested for LO events.  Include "
			    << "'set Factory:ShowerApproximation MEMatching' in input file."
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

        if ( !isMCatNLOHEvent )
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
	  SpinOneLorentzRotation rot =
	    intrinsicPtGenerator->kick(eventRecord().incoming(),
				       eventRecord().intermediates());
	  eventRecord().transform(rot);
	}
      }

      didRealign = realign();
    
      constituentReshuffle();

      break;

    } catch (RedoShower&) {

      resetWeights();

      if ( ++nTries > maxtry() )
	throw ShowerTriesVeto(maxtry());

      eventRecord().clear();
      eventRecord().prepare(sub,dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr()),pdfs());

      continue;

    } catch (...) {
      throw;
    }

  }

  return eventRecord().fillEventRecord(newStep(),firstInteraction(),didRealign);

}

void DipoleShowerHandler::constituentReshuffle() {

  if ( constituentReshuffler ) {
    constituentReshuffler->reshuffle(eventRecord().outgoing(),
				     eventRecord().incoming(),
				     eventRecord().intermediates());
  }

}

void DipoleShowerHandler::hardScales(Energy2 muf) {

  maxPt = generator()->maximumCMEnergy();

  if ( restrictPhasespace() ) {
    if ( !hardScaleIsMuF() || !firstInteraction() ) {
      if ( !eventRecord().outgoing().empty() ) {
	for ( PList::const_iterator p = eventRecord().outgoing().begin();
	      p != eventRecord().outgoing().end(); ++p )
	  maxPt = min(maxPt,(**p).momentum().mt());
      } else {
	assert(!eventRecord().hard().empty());
	Lorentz5Momentum phard(ZERO,ZERO,ZERO,ZERO);
	for ( PList::const_iterator p = eventRecord().hard().begin();
	      p != eventRecord().hard().end(); ++p )
	  phard += (**p).momentum();
	Energy mhard = phard.m();
	maxPt = mhard;
      }
      maxPt *= hardScaleFactor();
    } else {
      maxPt = hardScaleFactor()*sqrt(muf);
    }
    muPt = maxPt;
  } else {
    muPt = hardScaleFactor()*sqrt(muf);
  }

  for ( list<DipoleChain>::iterator ch = eventRecord().chains().begin();
	ch != eventRecord().chains().end(); ++ch ) {

    Energy minVetoScale = -1.*GeV;

    for ( list<Dipole>::iterator dip = ch->dipoles().begin();
	  dip != ch->dipoles().end(); ++dip ) {

      // max scale per config
      Energy maxFirst = 0.0*GeV;
      Energy maxSecond = 0.0*GeV;

      for ( vector<Ptr<DipoleSplittingKernel>::ptr>::iterator k =
	      kernels.begin(); k != kernels.end(); ++k ) {
	
	pair<bool,bool> conf = make_pair(true,false);

	if ( (**k).canHandle(dip->index(conf)) ) {
	  Energy scale =
	    evolutionOrdering()->hardScale(dip->emitter(conf),dip->spectator(conf),
					 dip->emitterX(conf),dip->spectatorX(conf),
					 **k,dip->index(conf));
	  maxFirst = max(maxFirst,scale);
	}

	conf = make_pair(false,true);

	if ( (**k).canHandle(dip->index(conf)) ) {
	  Energy scale =
	    evolutionOrdering()->hardScale(dip->emitter(conf),dip->spectator(conf),
					 dip->emitterX(conf),dip->spectatorX(conf),
					 **k,dip->index(conf));
	  maxSecond = max(maxSecond,scale);
	}

      }

      if ( dip->leftParticle()->vetoScale() >= ZERO ) {
	maxFirst = min(maxFirst,sqrt(dip->leftParticle()->vetoScale()));
	if ( minVetoScale >= ZERO )
	  minVetoScale = min(minVetoScale,sqrt(dip->leftParticle()->vetoScale()));
	else
	  minVetoScale = sqrt(dip->leftParticle()->vetoScale());
      }

      if ( dip->rightParticle()->vetoScale() >= ZERO ) {
	maxSecond = min(maxSecond,sqrt(dip->rightParticle()->vetoScale()));
	if ( minVetoScale >= ZERO )
	  minVetoScale = min(minVetoScale,sqrt(dip->rightParticle()->vetoScale()));
	else
	  minVetoScale = sqrt(dip->rightParticle()->vetoScale());
      }

      maxFirst = min(maxPt,maxFirst);
      dip->emitterScale(make_pair(true,false),maxFirst);

      maxSecond = min(maxPt,maxSecond);
      dip->emitterScale(make_pair(false,true),maxSecond);

    }

    if ( !evolutionOrdering()->independentDipoles() &&
	 chainOrderVetoScales &&
	 minVetoScale >= ZERO ) {
      for ( list<Dipole>::iterator dip = ch->dipoles().begin();
	    dip != ch->dipoles().end(); ++dip ) {
	dip->leftScale(min(dip->leftScale(),minVetoScale));
	dip->rightScale(min(dip->rightScale(),minVetoScale));
      }
    }

  }

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
       !doFSR ) {
    winner.didStopEvolving();
    return 0.0*GeV;
  }

  if ( index.initialStateEmitter() &&
       !doISR ) {
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

  pair<GeneratorMap::iterator,GeneratorMap::iterator> gens
    = generators().equal_range(candidate.index());

  Energy winnerScale = 0.0*GeV;
  GeneratorMap::iterator winnerGen = generators().end();

  for ( GeneratorMap::iterator gen = gens.first; gen != gens.second; ++gen ) {

    // (*) see NOTE above
    if ( !(gen->first == candidate.index()) )
      continue;

    if ( startScale <= gen->second->splittingKinematics()->IRCutoff() )
      continue;

    Energy dScale =
      gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),
						      spectator->momentum());

    // in very exceptional cases happening in DIS
    if ( isnan(dScale/GeV ) )
      throw RedoShower();

    candidate.scale(dScale);
    candidate.continuesEvolving();
    Energy hardScale = evolutionOrdering()->maxPt(startScale,candidate,*(gen->second->splittingKernel()));

    Energy maxPossible = 
      gen->second->splittingKinematics()->ptMax(candidate.scale(),
						candidate.emitterX(), candidate.spectatorX(),
						candidate.index(),
						*gen->second->splittingKernel());

    Energy ircutoff =
      optCutoff < gen->second->splittingKinematics()->IRCutoff() ?
      gen->second->splittingKinematics()->IRCutoff() :
      optCutoff;

    if ( maxPossible <= ircutoff ) {
      continue;
    }

    if ( maxPossible >= hardScale )
      candidate.hardPt(hardScale);
    else {
      hardScale = maxPossible;
      candidate.hardPt(maxPossible);
    }

    gen->second->generate(candidate,currentWeights_,optHardPt,optCutoff);
    Energy nextScale = evolutionOrdering()->evolutionScale(gen->second->lastSplitting(),*(gen->second->splittingKernel()));

    if ( nextScale > winnerScale ) {
      winner.fill(candidate);
      gen->second->completeSplitting(winner);
      winnerGen = gen;
      winnerScale = nextScale;
    }

    reweight_ *= gen->second->splittingWeight();

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
				    Energy optCutoff) {

  if ( nEmissions )
    if ( emDone == nEmissions )
      return;

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
      
      nextLeftScale = getWinner(dipoleWinner,*dip,make_pair(true,false),optHardPt,optCutoff);

      if ( nextLeftScale > winnerScale ) {
	winnerScale = nextLeftScale;
	winner = dipoleWinner;
	winnerDip = dip;
      }

      nextRightScale = getWinner(dipoleWinner,*dip,make_pair(false,true),optHardPt,optCutoff);

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
	  reweight_ *= w;
	}
      continue;
    }

    // otherwise perform the splitting

    didRadiate = true;

    isMCatNLOSEvent = false;
    isMCatNLOHEvent = false;

    pair<list<Dipole>::iterator,list<Dipole>::iterator> children;

    DipoleChain* firstChain = 0;
    DipoleChain* secondChain = 0;

    eventRecord().split(winnerDip,winner,children,firstChain,secondChain);

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
	reweight_ *= w;
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

    SpinOneLorentzRotation transform((inMomenta.first+inMomenta.second).findBoostToCM());

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

  for ( vector<Ptr<DipoleSplittingKernel>::ptr>::iterator k = kernels.begin();
	k != kernels.end(); ++k ) {
    if ( !(**k).alphaS() )
      (**k).alphaS(as);
    (**k).renormalizationScaleFreeze(theRenormalizationScaleFreeze);
    (**k).factorizationScaleFreeze(theFactorizationScaleFreeze);
  }

  // clear the generators to be rebuild
  // actually, there shouldn't be any generators
  // when this happens.
  generators().clear();

}

void DipoleShowerHandler::resetReweight(Ptr<DipoleSplittingReweight>::tptr rw) {
  for ( GeneratorMap::iterator k = generators().begin();
	k != generators().end(); ++k )
    k->second->splittingReweight(rw);
}

void DipoleShowerHandler::getGenerators(const DipoleIndex& ind,
					Ptr<DipoleSplittingReweight>::tptr rw) {

  bool gotone = false;

  for ( vector<Ptr<DipoleSplittingKernel>::ptr>::iterator k =
	  kernels.begin(); k != kernels.end(); ++k ) {
    if ( (**k).canHandle(ind) ) {

      if ( verbosity > 0 ) {
	generator()->log() << "DipoleShowerHandler encountered the dipole configuration\n"
			   << ind << " in event number "
			   << eventHandler()->currentEvent()->number()
			   << "\nwhich can be handled by the splitting kernel '"
			   << (**k).name() << "'.\n" << flush;
      }

      gotone = true;

      Ptr<DipoleSplittingGenerator>::ptr nGenerator =
	new_ptr(DipoleSplittingGenerator());
      nGenerator->doCompensate(theDoCompensate);
      nGenerator->splittingKernel(*k);
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
	  if ( (**k).canHandleEquivalent(ind,*(eq->second->splittingKernel()),eq->first) ) {

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

      generators().insert(make_pair(ind,nGenerator));


    }
  }

  if ( !gotone ) {
    generator()->logWarning(Exception() 
			    << "DipoleShowerHandler could not "
			    << "find a splitting kernel which is able "
			    << "to handle splittings off the dipole "
			    << ind << ".\n"
			    << "Please check the input files."
			    << Exception::warning);
  }

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void DipoleShowerHandler::doinit() {
  ShowerHandler::doinit();
  if ( theGlobalAlphaS )
    resetAlphaS(theGlobalAlphaS);
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
     << nEmissions << discardNoEmissions << firstMCatNLOEmission << doFSR << doISR
     << realignmentScheme << verbosity << printEvent
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV)
     << isMCatNLOSEvent << isMCatNLOHEvent << theShowerApproximation
     << theDoCompensate << theFreezeGrid << theDetuning
     << theEventReweight << theSplittingReweight << ounit(maxPt,GeV)
     << ounit(muPt,GeV);
}

void DipoleShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> kernels >> theEvolutionOrdering 
     >> constituentReshuffler >> intrinsicPtGenerator
     >> theGlobalAlphaS >> chainOrderVetoScales
     >> nEmissions >> discardNoEmissions >> firstMCatNLOEmission >> doFSR >> doISR
     >> realignmentScheme >> verbosity >> printEvent
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV)
     >> isMCatNLOSEvent >> isMCatNLOHEvent >> theShowerApproximation
     >> theDoCompensate >> theFreezeGrid >> theDetuning
     >> theEventReweight >> theSplittingReweight >> iunit(maxPt,GeV)
     >> iunit(muPt,GeV);
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


  static Switch<DipoleShowerHandler,bool> interfaceDoFSR
    ("DoFSR",
     "Switch on or off final state radiation.",
     &DipoleShowerHandler::doFSR, true, false, false);
  static SwitchOption interfaceDoFSROn
    (interfaceDoFSR,
     "On",
     "Switch on final state radiation.",
     true);
  static SwitchOption interfaceDoFSROff
    (interfaceDoFSR,
     "Off",
     "Switch off final state radiation.",
     false);

  static Switch<DipoleShowerHandler,bool> interfaceDoISR
    ("DoISR",
     "Switch on or off initial state radiation.",
     &DipoleShowerHandler::doISR, true, false, false);
  static SwitchOption interfaceDoISROn
    (interfaceDoISR,
     "On",
     "Switch on initial state radiation.",
     true);
  static SwitchOption interfaceDoISROff
    (interfaceDoISR,
     "Off",
     "Switch off initial state radiation.",
     false);

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
     "[experimental] Switch on or off the chain ordering for veto scales.",
     &DipoleShowerHandler::chainOrderVetoScales, true, false, false);
  static SwitchOption interfaceChainOrderVetoScalesOn
    (interfaceChainOrderVetoScales,
     "On",
     "Switch on chain ordering for veto scales.",
     true);
  static SwitchOption interfaceChainOrderVetoScalesOff
    (interfaceChainOrderVetoScales,
     "Off",
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
  static SwitchOption interfaceDiscardNoEmissionsOn
    (interfaceDiscardNoEmissions,
     "On",
     "Discard events without radiation.",
     true);
  static SwitchOption interfaceDiscardNoEmissionsOff 
    (interfaceDiscardNoEmissions,
     "Off",
     "Do not discard events without radiation.",
     false);

  interfaceDiscardNoEmissions.rank(-1);

  static Switch<DipoleShowerHandler,bool> interfaceFirstMCatNLOEmission
    ("FirstMCatNLOEmission",
     "[debug option] Only perform the first MC@NLO emission.",
     &DipoleShowerHandler::firstMCatNLOEmission, false, false, false);
  static SwitchOption interfaceFirstMCatNLOEmissionOn
    (interfaceFirstMCatNLOEmission,
     "On",
     "",
     true);
  static SwitchOption interfaceFirstMCatNLOEmissionOff 
    (interfaceFirstMCatNLOEmission,
     "Off",
     "",
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

}

