// -*- C++ -*-
//
// DipoleSplittingGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingGenerator class.
//
#include <config.h>
#include "DipoleSplittingGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/DipoleShower/DipoleShowerHandler.h"

using namespace Herwig;

DipoleSplittingGenerator::DipoleSplittingGenerator() 
  : HandlerBase(),
    theExponentialGenerator(0), prepared(false), presampling(false),
    theDoCompensate(false), theSplittingWeight(1.) {
  if ( ShowerHandler::currentHandler() )
    setGenerator(ShowerHandler::currentHandler()->generator());
}

DipoleSplittingGenerator::~DipoleSplittingGenerator() {
  if ( theExponentialGenerator ) {
    delete theExponentialGenerator;
    theExponentialGenerator = 0;
  }
}

IBPtr DipoleSplittingGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleSplittingGenerator::fullclone() const {
  return new_ptr(*this);
}

void DipoleSplittingGenerator::wrap(Ptr<DipoleSplittingGenerator>::ptr other) {
  assert(!prepared);
  theOtherGenerator = other;
}

void DipoleSplittingGenerator::resetVariations() {
  for ( map<string,double>::iterator w = currentWeights.begin();
	w != currentWeights.end(); ++w )
    w->second = 1.;
}

void DipoleSplittingGenerator::veto(const vector<double>&, double p, double r) {
  double factor = 1.;
  if ( splittingReweight() ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() && splittingReweight()->firstInteraction() ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() && splittingReweight()->secondaryInteractions() ) ) {
      factor = splittingReweight()->evaluate(generatedSplitting);
      theSplittingWeight *= (r-factor*p)/(r-p);
    }
  }
  splittingKernel()->veto(generatedSplitting, factor*p, r, currentWeights);
}

void DipoleSplittingGenerator::accept(const vector<double>&, double p, double r) {
  double factor = 1.;
  if ( splittingReweight() ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() && splittingReweight()->firstInteraction() ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() && splittingReweight()->secondaryInteractions() ) ) {
      factor = splittingReweight()->evaluate(generatedSplitting);
      theSplittingWeight *= factor;
    }
  }
  splittingKernel()->accept(generatedSplitting, factor*p, r, currentWeights);
}

void DipoleSplittingGenerator::prepare(const DipoleSplittingInfo& sp) {

  generatedSplitting = sp;

  generatedSplitting.splittingKinematics(splittingKernel()->splittingKinematics());
  generatedSplitting.splittingParameters().resize(splittingKernel()->nDimAdditional());

  if ( wrapping() ) {
    generatedSplitting.emitterData(theSplittingKernel->emitter(generatedSplitting.index()));  
    generatedSplitting.spectatorData(theSplittingKernel->spectator(generatedSplitting.index()));  
    generatedSplitting.emissionData(theSplittingKernel->emission(generatedSplitting.index()));  
    parameters.resize(theOtherGenerator->nDim());
    prepared = true;
    return;
  }

  generatedSplitting.emitterData(splittingKernel()->emitter(generatedSplitting.index()));  
  generatedSplitting.spectatorData(splittingKernel()->spectator(generatedSplitting.index()));  
  generatedSplitting.emissionData(splittingKernel()->emission(generatedSplitting.index()));  

  presampledSplitting = generatedSplitting;

  prepared = true;

  parameters.resize(nDim());

  theExponentialGenerator = 
    new exsample::exponential_generator<DipoleSplittingGenerator,UseRandom>();

  theExponentialGenerator->sampling_parameters().maxtry = maxtry();
  theExponentialGenerator->sampling_parameters().presampling_points = presamplingPoints();
  theExponentialGenerator->sampling_parameters().freeze_grid = freezeGrid();
  theExponentialGenerator->detuning(detuning());

  theExponentialGenerator->docompensate(theDoCompensate);
  theExponentialGenerator->function(this);
  theExponentialGenerator->initialize();

}

void DipoleSplittingGenerator::fixParameters(const DipoleSplittingInfo& sp,
					     Energy optHardPt) {

  assert(generator());

  assert(!presampling);
  assert(prepared);

  assert(sp.index() == generatedSplitting.index());

  generatedSplitting.scale(sp.scale());
  parameters[3] = sp.scale()/generator()->maximumCMEnergy();

  generatedSplitting.hardPt(sp.hardPt());

  parameters[0] = splittingKinematics()->ptToRandom(optHardPt == ZERO ? 
						    generatedSplitting.hardPt() : 
						    min(generatedSplitting.hardPt(),optHardPt),
						    sp.scale(),
						    sp.emitterX(), sp.spectatorX(),
						    generatedSplitting.index(),
						    *splittingKernel());

  size_t shift = 4;

  if ( generatedSplitting.index().emitterPDF().pdf() &&
       generatedSplitting.index().spectatorPDF().pdf() ) {
    generatedSplitting.emitterX(sp.emitterX());
    generatedSplitting.spectatorX(sp.spectatorX());
    parameters[4] = sp.emitterX();
    parameters[5] = sp.spectatorX();
    shift += 2;
  }

  if ( generatedSplitting.index().emitterPDF().pdf() &&
       !generatedSplitting.index().spectatorPDF().pdf() ) {
    generatedSplitting.emitterX(sp.emitterX());
    parameters[4] = sp.emitterX();
    ++shift;
  }

  if ( !generatedSplitting.index().emitterPDF().pdf() &&
       generatedSplitting.index().spectatorPDF().pdf() ) {
    generatedSplitting.spectatorX(sp.spectatorX());
    parameters[4] = sp.spectatorX();
    ++shift;
  }

  if ( splittingKernel()->nDimAdditional() )
    copy(sp.lastSplittingParameters().begin(),sp.lastSplittingParameters().end(),parameters.begin()+shift);

  if ( sp.emitter() )
    generatedSplitting.emitter(sp.emitter());

  if ( sp.spectator() )
    generatedSplitting.spectator(sp.spectator());

}

int DipoleSplittingGenerator::nDim() const {

  assert(!wrapping());
  assert(prepared);

  int ret = 4; // 0 pt, 1 z, 2 phi, 3 scale, 4/5 xs + parameters

  if ( generatedSplitting.index().emitterPDF().pdf() ) {
    ++ret;
  }  

  if ( generatedSplitting.index().spectatorPDF().pdf() ) {
    ++ret;
  }  

  ret += splittingKernel()->nDimAdditional();

  return ret;

}

const vector<bool>& DipoleSplittingGenerator::sampleFlags() {

  assert(!wrapping());

  if ( !theFlags.empty() )
    return theFlags;

  theFlags.resize(nDim(),false);
  theFlags[0] = true; theFlags[1] = true; theFlags[2] = true; // 0 pt, 1 z, 2 phi
  return theFlags;
}

const pair<vector<double>,vector<double> >& DipoleSplittingGenerator::support() {

  assert(!wrapping());

  if ( !theSupport.first.empty() )
    return theSupport;

  vector<double> lower(nDim(),0.);
  vector<double> upper(nDim(),1.);

  pair<double,double> kSupport = 
    generatedSplitting.splittingKinematics()->kappaSupport(generatedSplitting);

  pair<double,double> xSupport = 
    generatedSplitting.splittingKinematics()->xiSupport(generatedSplitting);

  lower[0] = kSupport.first;
  lower[1] = xSupport.first;

  upper[0] = kSupport.second;
  upper[1] = xSupport.second;

  theSupport.first = lower;
  theSupport.second = upper;

  return theSupport;

}

void DipoleSplittingGenerator::startPresampling() {
  assert(!wrapping());
  splittingKernel()->startPresampling(generatedSplitting.index());
  presampling = true;
}

void DipoleSplittingGenerator::stopPresampling() {
  assert(!wrapping());
  splittingKernel()->stopPresampling(generatedSplitting.index());
  presampling = false;
}

bool DipoleSplittingGenerator::haveOverestimate() const {

  assert(!wrapping());
  assert(prepared);

  return 
    generatedSplitting.splittingKinematics()->haveOverestimate() &&
    splittingKernel()->haveOverestimate(generatedSplitting);

}

bool DipoleSplittingGenerator::overestimate(const vector<double>& point) {

  assert(!wrapping());
  assert(prepared);
  assert(!presampling);
  assert(haveOverestimate());

  if ( ! generatedSplitting.splittingKinematics()->generateSplitting(point[0],point[1],point[2],
								     generatedSplitting,
								     *splittingKernel()) )
    return 0.;

  generatedSplitting.splittingKinematics()->prepareSplitting(generatedSplitting);

  return 
    ( generatedSplitting.splittingKinematics()->jacobianOverestimate() * 
      splittingKernel()->overestimate(generatedSplitting) );

}

double DipoleSplittingGenerator::invertOverestimateIntegral(double value) const {

  assert(!wrapping());
  assert(prepared);
  assert(!presampling);
  assert(haveOverestimate());

  return 
    splittingKernel()->invertOverestimateIntegral(generatedSplitting,value);

}

double DipoleSplittingGenerator::evaluate(const vector<double>& point) {

  assert(!wrapping());
  assert(prepared);
  assert(generator());

  DipoleSplittingInfo& split =
    ( !presampling ? generatedSplitting : presampledSplitting );

  split.continuesEvolving();

  size_t shift = 4;

  if ( presampling ) {

    split.scale(point[3] * generator()->maximumCMEnergy());

    if ( split.index().emitterPDF().pdf() &&
	 split.index().spectatorPDF().pdf() ) {
      split.emitterX(point[4]);
      split.spectatorX(point[5]);
      shift += 2;
    }

    if ( split.index().emitterPDF().pdf() &&
	 !split.index().spectatorPDF().pdf() ) {
      split.emitterX(point[4]);
      ++shift;
    }

    if ( !split.index().emitterPDF().pdf() &&
	 split.index().spectatorPDF().pdf() ) {
      split.spectatorX(point[4]);
      ++shift;
    }

    if ( splittingKernel()->nDimAdditional() )
      copy(point.begin()+shift,point.end(),split.splittingParameters().begin());

    split.hardPt(split.splittingKinematics()->ptMax(split.scale(),
						    split.emitterX(),
						    split.spectatorX(),
						    split.index(),
						    *splittingKernel()));

  }

  if ( ! split.splittingKinematics()->generateSplitting(point[0],point[1],point[2],split,*splittingKernel()) ) {
    split.lastValue(0.);
    return 0.;
  }

  split.splittingKinematics()->prepareSplitting(split);

  if ( split.stoppedEvolving() ) {
    split.lastValue(0.);
    return 0.;
  }

  if ( !presampling )
    splittingKernel()->clearAlphaPDFCache();
  double kernel = splittingKernel()->evaluate(split);
  double jac = split.splittingKinematics()->jacobian();

  // multiply in the profile scales when relevant
  assert(ShowerHandler::currentHandler());
  if ( ShowerHandler::currentHandler()->firstInteraction() &&
       ShowerHandler::currentHandler()->profileScales() &&
       !presampling ) {
    Energy hard = ShowerHandler::currentHandler()->hardScale();
    if ( hard > ZERO )
      kernel *= ShowerHandler::currentHandler()->profileScales()->hardScaleProfile(hard,split.lastPt());
  }

  split.lastValue( abs(jac) * kernel );

  if ( isnan(split.lastValue()) || isinf(split.lastValue()) ) {
    generator()->log() << "DipoleSplittingGenerator:evaluate(): problematic splitting kernel encountered for "
		       << splittingKernel()->name() << "\n" << flush;
    split.lastValue(0.0);
  }

  if ( kernel < 0. )
    return 0.;

  return split.lastValue();

}

void DipoleSplittingGenerator::doGenerate(map<string,double>& variations,
					  Energy optCutoff) {

  assert(!wrapping());

  double res = 0.;

  Energy startPt = generatedSplitting.hardPt();
  double optKappaCutoff = 0.0;
  if ( optCutoff > splittingKinematics()->IRCutoff() ) {
    optKappaCutoff = splittingKinematics()->ptToRandom(optCutoff,
						       generatedSplitting.scale(),
						       generatedSplitting.emitterX(), 
						       generatedSplitting.spectatorX(),
						       generatedSplitting.index(),
						       *splittingKernel());
  }

  resetVariations();
  theSplittingWeight = 1.;

  while (true) {
    try {
      if ( optKappaCutoff == 0.0 ) {
	res = theExponentialGenerator->generate();
      } else {
	res = theExponentialGenerator->generate(optKappaCutoff);
      }
    } catch (exsample::exponential_regenerate&) {
      resetVariations();
      theSplittingWeight = 1.;
      generatedSplitting.hardPt(startPt);
      continue;
    } catch (exsample::hit_and_miss_maxtry&) {
      throw DipoleShowerHandler::RedoShower();
    } catch (exsample::selection_maxtry&) {
      throw DipoleShowerHandler::RedoShower();
    }
    break;
  }

  for ( map<string,double>::const_iterator w = currentWeights.begin();
	w != currentWeights.end(); ++w ) {
    map<string,double>::iterator v = variations.find(w->first);
    if ( v != variations.end() )
      v->second *= w->second;
    else
      variations[w->first] = w->second;
  }

  if ( res == 0. ) {
    generatedSplitting.lastPt(0.0*GeV);
    generatedSplitting.didStopEvolving();
  } else {

    generatedSplitting.continuesEvolving();

    if ( theMCCheck )
      theMCCheck->book(generatedSplitting.emitterX(),
		       generatedSplitting.spectatorX(),
		       generatedSplitting.scale(),
		       startPt,
		       generatedSplitting.lastPt(),
		       generatedSplitting.lastZ(),
		       1.);

  }

}

Energy DipoleSplittingGenerator::generate(const DipoleSplittingInfo& split,
					  map<string,double>& variations,
					  Energy optHardPt,
					  Energy optCutoff) {

  fixParameters(split,optHardPt);

  if ( wrapping() ) {
    return theOtherGenerator->generateWrapped(generatedSplitting,variations,optHardPt,optCutoff);
  }

  doGenerate(variations,optCutoff);

  return generatedSplitting.lastPt();

}

Energy DipoleSplittingGenerator::generateWrapped(DipoleSplittingInfo& split,
						 map<string,double>& variations,
						 Energy optHardPt,
						 Energy optCutoff) {

  assert(!wrapping());

  DipoleSplittingInfo backup = generatedSplitting;
  generatedSplitting = split;

  fixParameters(split,optHardPt);

  try {
    doGenerate(variations,optCutoff);
  } catch (...) {
    split = generatedSplitting;
    generatedSplitting = backup;
    throw;
  }

  Energy pt = generatedSplitting.lastPt();

  split = generatedSplitting;
  generatedSplitting = backup;

  return pt;

}

void DipoleSplittingGenerator::completeSplitting(DipoleSplittingInfo& sp) const {
  pair<bool,bool> conf = sp.configuration();
  sp = generatedSplitting;
  sp.configuration(conf);
}

Ptr<DipoleSplittingKernel>::tptr DipoleSplittingGenerator::splittingKernel() const { 
  if ( wrapping() )
    return theOtherGenerator->splittingKernel();
  return theSplittingKernel;
}

Ptr<DipoleSplittingReweight>::tptr DipoleSplittingGenerator::splittingReweight() const { 
  if ( wrapping() )
    return theOtherGenerator->splittingReweight();
  return theSplittingReweight;
}

Ptr<DipoleSplittingKinematics>::tptr DipoleSplittingGenerator::splittingKinematics() const { 
  if ( wrapping() )
    return theOtherGenerator->splittingKinematics();
  return theSplittingKernel->splittingKinematics();
}

void DipoleSplittingGenerator::splittingKernel(Ptr<DipoleSplittingKernel>::tptr sp) { 
  theSplittingKernel = sp;
  if ( theSplittingKernel->mcCheck() )
    theMCCheck = theSplittingKernel->mcCheck();
}

void DipoleSplittingGenerator::splittingReweight(Ptr<DipoleSplittingReweight>::tptr sp) { 
  theSplittingReweight = sp;
}

void DipoleSplittingGenerator::debugGenerator(ostream& os) const {

  os << "--- DipoleSplittingGenerator ---------------------------------------------------\n";

  os << " generating splittings using\n"
     << " splittingKernel = " << splittingKernel()->name()
     << " splittingKinematics = " << generatedSplitting.splittingKinematics()->name() << "\n"
     << " to sample splittings of type:\n";

  os << generatedSplitting;

  os << "--------------------------------------------------------------------------------\n";

}

void DipoleSplittingGenerator::debugLastEvent(ostream& os) const {

  os << "--- DipoleSplittingGenerator ---------------------------------------------------\n";

  os << " last generated event:\n";

  os << generatedSplitting;

  os << "--------------------------------------------------------------------------------\n";

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << theOtherGenerator << theSplittingKernel << theSplittingReweight << theMCCheck << theDoCompensate;
}

void DipoleSplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theOtherGenerator >> theSplittingKernel >> theSplittingReweight >> theMCCheck >> theDoCompensate;
}

ClassDescription<DipoleSplittingGenerator> DipoleSplittingGenerator::initDipoleSplittingGenerator;
// Definition of the static class description member.

void DipoleSplittingGenerator::Init() {

  static ClassDocumentation<DipoleSplittingGenerator> documentation
    ("DipoleSplittingGenerator is used by the dipole shower "
     "to sample splittings from a given dipole splitting kernel.");


  static Reference<DipoleSplittingGenerator,DipoleSplittingKernel> interfaceSplittingKernel
    ("SplittingKernel",
     "Set the splitting kernel to sample from.",
     &DipoleSplittingGenerator::theSplittingKernel, false, false, true, false, false);

  static Reference<DipoleSplittingGenerator,DipoleSplittingReweight> interfaceSplittingReweight
    ("SplittingReweight",
     "Set the splitting reweight.",
     &DipoleSplittingGenerator::theSplittingReweight, false, false, true, true, false);

  static Reference<DipoleSplittingGenerator,DipoleMCCheck> interfaceMCCheck
    ("MCCheck",
     "[debug option] MCCheck",
     &DipoleSplittingGenerator::theMCCheck, false, false, true, true, false);

  interfaceMCCheck.rank(-1);

}

