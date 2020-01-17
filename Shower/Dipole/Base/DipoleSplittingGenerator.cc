// -*- C++ -*-
//
// DipoleSplittingGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"
#include "ThePEG/Repository/UseRandom.h"

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
    if ( ( ShowerHandler::currentHandler()->firstInteraction() &&
          splittingReweight()->firstInteraction() ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() &&
	   splittingReweight()->secondaryInteractions() ) ) {
      if ( !splittingReweight()->hintOnly(generatedSplitting) ) {
	factor = 
	  splittingReweight()->evaluate(generatedSplitting)/splittingReweight()->hint(generatedSplitting);
	theSplittingWeight *= (r-factor*p)/(r-p);
	theSplittingWeightVector.push_back(std::make_tuple(generatedSplitting.lastPt(),(r-factor*p)/(r-p),false));
      }
    }
  }
  splittingKernel()->veto(generatedSplitting, factor*p, r, currentWeights);
}

void DipoleSplittingGenerator::accept(const vector<double>&, double p, double r) {
  double factor = 1.;
  if ( splittingReweight() ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() &&
          splittingReweight()->firstInteraction() ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() &&
	   splittingReweight()->secondaryInteractions() ) ) {
      if ( !splittingReweight()->hintOnly(generatedSplitting) ) {
	factor = 
	  splittingReweight()->evaluate(generatedSplitting)/splittingReweight()->hint(generatedSplitting);
	theSplittingWeight *= factor;
        theSplittingWeightVector.push_back(std::make_tuple(generatedSplitting.lastPt(),factor,true));
      } else {
	theSplittingWeightVector.push_back(std::make_tuple(generatedSplitting.lastPt(),1.0,true));
      }
    }
  }
  splittingKernel()->accept(generatedSplitting, factor*p, r, currentWeights);
}

void DipoleSplittingGenerator::prepare(const DipoleSplittingInfo& sp) {

  generatedSplitting = sp;

  generatedSplitting.splittingKinematics(splittingKernel()->splittingKinematics());
  // The splitting kernel is needed for spin correlations
  generatedSplitting.splittingKernel(splittingKernel());
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

  // If dealing with a decay, need to set recoilMass
  if ( generatedSplitting.index().incomingDecaySpectator() ||    
       generatedSplitting.index().incomingDecayEmitter() ) 
    generatedSplitting.recoilMass(sp.recoilMass());

  // Need to copy emitter and spectator masses
  generatedSplitting.emitterMass(sp.emitterMass());
  generatedSplitting.spectatorMass(sp.spectatorMass());
  
  // Counter to track if there is an off-shell
  // emitter AND/OR spectator
  int count = parameters.size()-1;    

  // Off shell spectator mass
  if ( sp.index().offShellSpectator() ) {
    parameters[count] = sp.spectatorMass()/generator()->maximumCMEnergy();
    count -= 1;
  }
  
  // Off shell emitter mass
  if ( sp.index().offShellEmitter() )
    parameters[count] = sp.emitterMass()/generator()->maximumCMEnergy();

  
  // If not a decay, point[3] samples over the dipole scale
  if ( !sp.index().incomingDecaySpectator() && !sp.index().incomingDecayEmitter() )
    parameters[3] = sp.scale()/generator()->maximumCMEnergy();
  // If it is a decay, point[3] samples over the recoilMass
  else 
    parameters[3] = sp.recoilMass()/generator()->maximumCMEnergy();
  
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
    copy(sp.lastSplittingParameters().begin(),
         sp.lastSplittingParameters().end(),
         parameters.begin()+shift);

  if ( sp.emitter() )
    generatedSplitting.emitter(sp.emitter());

  if ( sp.spectator() )
    generatedSplitting.spectator(sp.spectator());

}

int DipoleSplittingGenerator::nDim() const {

  assert(!wrapping());
  assert(prepared);

  // Note this use of [3] for either the scale or the recoil mass
  // is a bit of a nasty hack.
  int ret = 4; // 0 pt, 1 z, 2 phi, 3 scale or recoilMass, 4/5 xs + parameters

  if ( generatedSplitting.index().emitterPDF().pdf() ) {
    ++ret;
  }  

  if ( generatedSplitting.index().spectatorPDF().pdf() ) {
    ++ret;
  }  

  ret += splittingKernel()->nDimAdditional();
  assert(splittingKernel()->nDimAdditional() == 0);

  // Put off-shell spectator mass at back [-1]
  // followed by off-shell emitter mass (i.e. [-1] or [-2])

  // Off-shell emitter
  if ( generatedSplitting.index().offShellEmitter() )
    ++ret;
  
  // Off-shell spectator
  if ( generatedSplitting.index().offShellSpectator() )
    ++ret;
  
  
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

double DipoleSplittingGenerator::overestimate(const vector<double>& point) {

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

    // Counter to track if there is an off-shell
    // emitter AND/OR spectator
    int count = parameters.size()-1;    

    // Sample over off-shell emitter and spectator masss
    // Do not sample if zero mass or off-shell
    if ( split.index().spectatorData()->mass() != ZERO ) {
      if ( !split.index().offShellSpectator() )
	split.spectatorMass(split.index().spectatorData()->mass());
      else {
	split.spectatorMass(point[count] * generator()->maximumCMEnergy());
	count -= 1;
      }
    }

    if ( split.index().emitterData()->mass() != ZERO ) {
      if ( !split.index().offShellEmitter() ) 
	split.emitterMass(split.index().emitterData()->mass());
      else
	split.emitterMass(point[count] * generator()->maximumCMEnergy());
    }
    
    // If not a decay, point[3] samples over the dipole scale
    if ( ! split.index().incomingDecaySpectator()
	 && ! split.index().incomingDecayEmitter() )
      split.scale(point[3] * generator()->maximumCMEnergy());
      
    // For dipoles containing a decayed spectator:
    // 1) Use point[3] to sample over the recoil mass
    // 2) The dipole scale is the spectator mass
    else if ( split.index().incomingDecaySpectator() ) {
      split.recoilMass(point[3] * generator()->maximumCMEnergy());
      assert(split.spectatorMass() != ZERO );
      split.scale(split.spectatorMass());
    }
    // Not currently intended to work with decaying emitters
    else
      assert(false);
    
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
						    split,
						    *splittingKernel()));
  }

  if ( ! split.splittingKinematics()->generateSplitting(point[0],
                                                        point[1],
                                                        point[2],
                                                        split,
                                                        *splittingKernel()) ) {
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
      kernel *= ShowerHandler::currentHandler()->profileScales()->
                hardScaleProfile(hard,split.lastPt());
  }

  split.lastValue( abs(jac) * kernel );

  if ( ! isfinite(split.lastValue()) ) {
    generator()->log() << "DipoleSplittingGenerator:evaluate():"
               <<"problematic splitting kernel encountered for "
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
  double enhance = 1.;
  bool detuningOff = false;
  if ( splittingReweight() ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() &&
          splittingReweight()->firstInteraction() ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() &&
	   splittingReweight()->secondaryInteractions() ) ) {
      enhance = splittingReweight()->hint(generatedSplitting);
      if ( splittingReweight()->hintOnly(generatedSplitting) )
	detuningOff = true;
    }
  }

  bool hintOnly = false;
  if ( splittingReweight() ) hintOnly = splittingReweight()->hintOnly(generatedSplitting);
  while (true) {
    theExponentialGenerator->detuning(detuning());
    if ( detuningOff )
      theExponentialGenerator->detuning(1.0);
    try {
      if ( optKappaCutoff == 0.0 ) {
        theSplittingWeightVector.clear();
	res = theExponentialGenerator->generate(enhance);
      } else {
	theSplittingWeightVector.clear();
	res = theExponentialGenerator->generate(optKappaCutoff,enhance);
      }
      //Partial unweighting
      if ( partialUnweighting && !hintOnly ) {
	if ( abs(theSplittingWeight)/theReferenceWeight < 1.0 ) {
	  double r = UseRandom::rnd(1.0);
	  if ( abs(theSplittingWeight)/theReferenceWeight < r ) {
	    theSplittingWeight = 1.;
	    continue;
	  } else {
	    theSplittingWeight = theSplittingWeight/abs(theSplittingWeight)*theReferenceWeight;
	  }
	}
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

double DipoleSplittingGenerator::sudakovExpansion(const DipoleSplittingInfo& split,
                                                  Energy down,Energy fixedScale){
  fixParameters(split);
  if ( wrapping() ) {
    return theOtherGenerator->wrappedSudakovExpansion( generatedSplitting, down,fixedScale);
  }
  return dosudakovExpansion( split, down,fixedScale);
}

double DipoleSplittingGenerator::sudakov(const DipoleSplittingInfo& split,Energy down){
  fixParameters(split);
  if ( wrapping() ) {
    return theOtherGenerator->wrappedSudakov( generatedSplitting, down);
  }
  return dosudakov( split, down);
}



double DipoleSplittingGenerator::dosudakovExpansion(const DipoleSplittingInfo& ,
                                                    Energy down,Energy fixedScale){
    assert(down > splittingKinematics()->IRCutoff());
  
  double optKappaCutoffd =
      splittingKinematics()->ptToRandom(down,
                                        generatedSplitting.scale(),
                                        generatedSplitting.emitterX(),
                                        generatedSplitting.spectatorX(),
                                        generatedSplitting.index(),
                                        *splittingKernel());
  
  
  double optKappaCutoffu = splittingKinematics()->ptToRandom(generatedSplitting.hardPt(),
                                                             generatedSplitting.scale(),
                                                             generatedSplitting.emitterX(),
                                                             generatedSplitting.spectatorX(),
                                                             generatedSplitting.index(),
                                                             *splittingKernel());
  
  
  pair<double,double> xSupport =
    generatedSplitting.splittingKinematics()->xiSupport(generatedSplitting);

  
    vector<double> RN;
    RN.resize(3);
    double res=0.;
    double resq=0.;
    double varx=10.;
    int k=0;
  
    generatedSplitting.setCalcFixedExpansion(true);
    generatedSplitting.fixedScale(fixedScale);
    
    while (  k<1000 ){
      k+=1.;
      RN[0]= optKappaCutoffd+(optKappaCutoffu-optKappaCutoffd)*UseRandom::rnd(); //PT
      RN[1]=xSupport.first+UseRandom::rnd()*(xSupport.second-xSupport.first); //
      RN[2]= UseRandom::rnd(); //PHI
      double tmp=(xSupport.second-xSupport.first)*
		 (optKappaCutoffu-optKappaCutoffd)*
                 evaluate(RN);
      res+= tmp;
      resq+=pow(tmp,2.);
      if(k%50==0.){
	varx=sqrt((resq/pow(1.*k,2)-pow(res,2)/pow(1.*k,3)))/(res/(1.0*k));
        if(varx<theSudakovAccuracy)break;
      }
    }
    generatedSplitting.setCalcFixedExpansion(false);

    return -res/(1.0*k);
}



double DipoleSplittingGenerator::dosudakov(const DipoleSplittingInfo& ,Energy down){
  
    

  double optKappaCutoffd = splittingKinematics()->ptToRandom(down,
                                                            generatedSplitting.scale(),
                                                            generatedSplitting.emitterX(),
                                                            generatedSplitting.spectatorX(),
                                                            generatedSplitting.index(),
                                                            *splittingKernel());
  
  
  double optKappaCutoffu = splittingKinematics()->ptToRandom(generatedSplitting.hardPt(),
                                                            generatedSplitting.scale(),
                                                            generatedSplitting.emitterX(),
                                                            generatedSplitting.spectatorX(),
                                                            generatedSplitting.index(),
                                                            *splittingKernel());
 

  
#ifndef NDEBUG
  pair<double,double> kSupport =
    generatedSplitting.splittingKinematics()->kappaSupport(generatedSplitting);
#endif
  assert(kSupport.first==0&&kSupport.second==1);

  pair<double,double> xSupport =
    generatedSplitting.splittingKinematics()->xiSupport(generatedSplitting);

 
    vector<double> RN;
    RN.resize(3);
  
    double res=0.;
    double resq=0.;
    double var=10.;
    double varx=10.;
    int k=0;
  while (((k<40.||var>theSudakovAccuracy)&&k<50000)){
    k+=1.;
    RN[0]= optKappaCutoffd+(optKappaCutoffu-optKappaCutoffd)*UseRandom::rnd(); //PT
    RN[1]=xSupport.first+UseRandom::rnd()*(xSupport.second-xSupport.first); //Z
    RN[2]=UseRandom::rnd(); //PHI
    double tmp=(xSupport.second-xSupport.first)*
	       (optKappaCutoffu-optKappaCutoffd)*
               evaluate(RN);
    
    res+= tmp;
    resq+=pow(tmp,2.);
    if(k%20==0.){
      varx=sqrt((resq/pow(1.*k,2)-pow(res,2)/pow(1.*k,3)));
      var=  (exp(-(res)/(1.0*k)+varx)-exp(-(res)/(1.0*k)-varx))/exp(-res/(1.0*k));
    }      

  }
 
  
    return exp(-res/(1.0*k));
    
}


double DipoleSplittingGenerator::wrappedSudakovExpansion(DipoleSplittingInfo& split,
                                                Energy down,Energy fixedScale) {
  
  assert(!wrapping());
  
  DipoleSplittingInfo backup = generatedSplitting;
  generatedSplitting = split;
  
  fixParameters(split);
  double res=dosudakovExpansion( split, down,fixedScale);
  
  split = generatedSplitting;
  generatedSplitting = backup;
  
  return res;
  
}

double DipoleSplittingGenerator::wrappedSudakov(DipoleSplittingInfo& split,
						 Energy down) {

  assert(!wrapping());

  DipoleSplittingInfo backup = generatedSplitting;
  generatedSplitting = split;
  
  fixParameters(split);
  double res=dosudakov( split, down);

  split = generatedSplitting;
  generatedSplitting = backup;

  return res;

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

