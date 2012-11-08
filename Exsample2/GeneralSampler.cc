// -*- C++ -*-
//
// GeneralSampler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralSampler class.
//

#include "GeneralSampler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/LoopGuard.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <boost/progress.hpp>

using namespace Herwig;

GeneralSampler::GeneralSampler() 
  : theVerbose(false), theFlatSubprocesses(false), 
    isSampling(false), theUpdateAfter(1),
    crossSectionCalls(0), gotCrossSections(false),
    theIntegratedXSec(0.), theIntegratedXSecErr(0.),
    theSumWeights(0.), norm(0.) {}

GeneralSampler::~GeneralSampler() {}

IBPtr GeneralSampler::clone() const {
  return new_ptr(*this);
}

IBPtr GeneralSampler::fullclone() const {
  return new_ptr(*this);
}

double sign(double x) {
  return x >= 0. ? 1. : -1.;
}

void GeneralSampler::initialize() {

  if ( theBinSampler->isUnweighting() && eventHandler()->weighted() ) {
    throw InitException() << "weighted events requested from unweighting bin sampler object.";
  }

  if ( !samplers.empty() )
    return;

  boost::progress_display* progressBar = 0;
  if ( !theVerbose ) {
    cout << "integrating subprocesses";
    progressBar = new boost::progress_display(eventHandler()->nBins(),cout);
  }

  for ( int b = 0; b < eventHandler()->nBins(); ++b ) {
    Ptr<BinSampler>::ptr s = theBinSampler->cloneMe();
    s->eventHandler(eventHandler());
    s->bin(b);
    lastSampler = s;
    s->initialize(theVerbose);
    samplers[b] = s;
    if ( !theVerbose )
      ++(*progressBar);
    if ( s->nanPoints() && theVerbose ) {
      cout << "warning: " 
	   << s->nanPoints() << " of "
	   << s->allPoints() << " points with nan or inf weight.\n"
	   << flush;
    }
  }

  updateCrossSections(true);

  if ( samplers.empty() ) {
    throw Exception() << "No processes with non-zero cross section present."
		      << Exception::abortnow;
  }

  if ( theVerbose )
    cout << "estimated total cross section is ( "
	 << integratedXSec()/nanobarn << " +/- "
	 << integratedXSecErr()/nanobarn << " ) nb\n" << flush;

  if ( progressBar )
    delete progressBar;

}

double GeneralSampler::generate() {

  long tries = 0;
  long excptTries = 0;

  gotCrossSections = false;

  if ( !theFlatSubprocesses )
    lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
  else {
    map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
    advance(s,(size_t)(UseRandom::rnd()*samplers.size()));
    lastSampler = s->second;
  }

  while ( true ) {

    try {
      lastSampler->generate(eventHandler()->weighted());
    } catch (BinSampler::NewMaximum& update) {
      if ( !eventHandler()->weighted() ) {
	unsigned long skip = 
	  (unsigned long)(lastSampler->acceptedPoints()*(1.-update.oldMaxWeight/update.newMaxWeight));
	map<Ptr<BinSampler>::tptr,unsigned long>::iterator s = skipMap.find(lastSampler);
	if ( s != skipMap.end() )
	  s->second += skip;
	else if ( skip != 0 )
	  skipMap[lastSampler] = skip;
	lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
	tries = 0;
	if ( ++excptTries == eventHandler()->maxLoop() )
	  break;
	continue;
      }
    } catch(BinSampler::UpdateCrossSections) {
      updateCrossSections();
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      tries = 0;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    } catch (...) {
      throw;
    }

    if ( isnan(lastSampler->lastWeight()) || isinf(lastSampler->lastWeight()) ) {
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      tries = 0;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    if ( eventHandler()->weighted() && lastSampler->lastWeight() == 0.0 ) {
      lastSampler->accept();
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      tries = 0;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    if ( eventHandler()->weighted() || lastSampler->isUnweighting() )
      break;

    if ( abs(lastSampler->lastWeight())/lastSampler->maxWeight() > UseRandom::rnd() ) {
      if ( skipMap.empty() )
	break;
      map<Ptr<BinSampler>::tptr,unsigned long>::iterator s = skipMap.find(lastSampler);
      if ( s == skipMap.end() )
	break;
      s->second -= 1;
      if ( s->second == 0 )
	skipMap.erase(s);
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      tries = 0;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    if ( ++tries == eventHandler()->maxLoop() ) {
      throw MaxTryException()
	<< "Maximum number of unweighting tries reached in GeneralSampler::generate()\n"
	<< "for process " << lastSampler->process()
	<< Exception::eventerror;
    }

  }

  if ( excptTries == eventHandler()->maxLoop() )
    throw Exception()
      << "GeneralSampler::generate() : Maximum number of tries to re-run event "
      << "selection reached. Aborting now." << Exception::runerror;

  lastPoint() = lastSampler->lastPoint();
  lastSampler->accept();

  if ( !eventHandler()->weighted() ) {
    theSumWeights += sign(lastSampler->lastWeight());
    return sign(lastSampler->lastWeight());
  } else {
    double w = lastSampler->lastWeight()/(norm*lastSampler->bias());
    theSumWeights += w;
    return w;
  }
  return 0.;

}

void GeneralSampler::rejectLast() {
  lastSampler->reject();
  if ( !eventHandler()->weighted() ) {
    theSumWeights -= sign(lastSampler->lastWeight());
  } else {
    theSumWeights -= 
      lastSampler->lastWeight()/(norm*lastSampler->bias());
  }
}

void GeneralSampler::currentCrossSections() const {

  if ( !isSampling )
    return;

  if ( gotCrossSections )
    return;

  if ( ++crossSectionCalls == theUpdateAfter ) {
    crossSectionCalls = 0;
  } else return;

  double xsec = 0.;
  double var = 0.;

  for ( map<double,Ptr<BinSampler>::ptr>::const_iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    xsec += s->second->averageWeight();
    var += s->second->averageWeightVariance();
  }

  theIntegratedXSec = xsec;
  theIntegratedXSecErr = sqrt(var);

  gotCrossSections = true;

}

void GeneralSampler::updateCrossSections(bool) {

  double xsec = 0.;
  double var = 0.;
  double sumbias = 0.;

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    if ( (isSampling && s->second == lastSampler) ||
	 !isSampling )
      s->second->nextIteration();
    if ( isSampling && s->second == lastSampler ) {
      s->second->maxWeight(s->second->iterations().back().maxWeight());
      s->second->minWeight(s->second->iterations().back().minWeight());
    }
    xsec += s->second->averageWeight();
    var += s->second->averageWeightVariance();
    sumbias += s->second->averageAbsWeight();
  }

  theIntegratedXSec = xsec;
  theIntegratedXSecErr = sqrt(var);
  norm = sumbias;

  if ( sumbias == 0.0 ) {
    samplers.clear();
    theIntegratedXSec = ZERO;
    theIntegratedXSecErr = ZERO;
    return;
  }

  map<double,Ptr<BinSampler>::ptr> newSamplers;
  double current = 0.;

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    double abssw = s->second->averageAbsWeight();
    if ( abssw == 0.0 )
      continue;
    s->second->bias(abssw/sumbias);
    current += abssw;
    newSamplers[current/sumbias] = s->second;
  }

  samplers = newSamplers;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void GeneralSampler::dofinish() {
  set<string> compensating;
  for ( map<double,Ptr<BinSampler>::ptr>::const_iterator s =
	  samplers.begin(); s != samplers.end(); ++s ) {
    if ( s->second->compensating() ||
	 skipMap.find(s->second) != skipMap.end() ) {
      compensating.insert(s->second->process());
    }
    if ( s->second->nanPoints() ) {
      generator()->logWarning(Exception()
			      << "warning: " 
			      << s->second->nanPoints() << " of "
			      << s->second->allPoints() << " points with nan or inf weight\n"
			      << "in " << s->second->process() << Exception::warning);
      if ( theVerbose ) {
	cout << "warning: " 
	     << s->second->nanPoints() << " of "
	     << s->second->allPoints() << " points with nan or inf weight\n"
	     << "in " << s->second->process() << "\n";
      }
    }
    s->second->finalize(theVerbose);
  }
  updateCrossSections();
  if ( theVerbose ) {
    if ( !compensating.empty() ) {
      cout << "warning: sampling for the following processes is still compensating:\n";
      for ( set<string>::const_iterator c = compensating.begin();
	    c != compensating.end(); ++c )
	cout << *c << "\n";
    }
    cout << "final integrated cross section is ( "
	 << integratedXSec()/nanobarn << " +/- "
	 << integratedXSecErr()/nanobarn << " ) nb\n" << flush;
  }
  if ( !compensating.empty() ) {
    generator()->logWarning(Exception()
			    << "Warning: Some samplers are still in compensating mode."
			    << Exception::warning);
  }
  SamplerBase::dofinish();
}

void GeneralSampler::doinitrun() {
  SamplerBase::doinitrun();
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    s->second->eventHandler(eventHandler());
    s->second->initialize(false);
    assert( !s->second->iterations().empty() );
    s->second->maxWeight(s->second->iterations().back().maxWeight());
    s->second->minWeight(s->second->iterations().back().minWeight());
  }
  isSampling = true;
}

void GeneralSampler::rebind(const TranslationMap & trans) {
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = 
	  samplers.begin(); s != samplers.end(); ++s )
    s->second = trans.translate(s->second);
  SamplerBase::rebind(trans);
}

IVector GeneralSampler::getReferences() {
  IVector ret = SamplerBase::getReferences();
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = 
	  samplers.begin(); s != samplers.end(); ++s )
    ret.push_back(s->second);
  return ret;
}

void GeneralSampler::persistentOutput(PersistentOStream & os) const {
  os << theBinSampler << theVerbose << theFlatSubprocesses 
     << samplers << lastSampler << theUpdateAfter
     << theIntegratedXSec << theIntegratedXSecErr
     << norm;
}

void GeneralSampler::persistentInput(PersistentIStream & is, int) {
  is >> theBinSampler >> theVerbose >> theFlatSubprocesses 
     >> samplers >> lastSampler >> theUpdateAfter
     >> theIntegratedXSec >> theIntegratedXSecErr
     >> norm;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<GeneralSampler,SamplerBase>
  describeHerwigGeneralSampler("Herwig::GeneralSampler", "HwExsample2.so");

void GeneralSampler::Init() {

  static ClassDocumentation<GeneralSampler> documentation
    ("A GeneralSampler class");


  static Reference<GeneralSampler,BinSampler> interfaceBinSampler
    ("BinSampler",
     "The bin sampler to be used.",
     &GeneralSampler::theBinSampler, false, false, true, false, false);


  static Parameter<GeneralSampler,size_t> interfaceUpdateAfter
    ("UpdateAfter",
     "Update cross sections every number of events.",
     &GeneralSampler::theUpdateAfter, 1, 1, 0,
     false, false, Interface::lowerlim);


  static Switch<GeneralSampler,bool> interfaceVerbose
    ("Verbose",
     "",
     &GeneralSampler::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseOn
    (interfaceVerbose,
     "On",
     "",
     true);
  static SwitchOption interfaceVerboseOff
    (interfaceVerbose,
     "Off",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceFlatSubprocesses
    ("FlatSubprocesses",
     "[debug] ",
     &GeneralSampler::theFlatSubprocesses, false, false, false);
  static SwitchOption interfaceFlatSubprocessesOn
    (interfaceFlatSubprocesses,
     "On",
     "",
     true);
  static SwitchOption interfaceFlatSubprocessesOff
    (interfaceFlatSubprocesses,
     "Off",
     "",
     false);

  interfaceFlatSubprocesses.rank(-1);

}

