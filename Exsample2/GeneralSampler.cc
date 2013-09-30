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

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <boost/progress.hpp>

using namespace Herwig;

GeneralSampler::GeneralSampler() 
  : theVerbose(false), theFlatSubprocesses(false), 
    isSampling(false), theUpdateAfter(1),
    crossSectionCalls(0), gotCrossSections(false),
    theIntegratedXSec(0.), theIntegratedXSecErr(0.),
    theSumWeights(0.), theSumWeights2(0.), 
    theAttempts(0), theAccepts(0),
    norm(0.), runCombinationData(false),
    theMaxWeight(0.0), theAlmostUnweighted(false) {}

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

  if ( !theBinSampler->isUnweighting() && !eventHandler()->weighted() ) {
    throw InitException() << "unweighted events requested from weighted bin sampler object.";
  }

  if ( theBinSampler->isUnweighting() && theFlatSubprocesses ) {
    throw InitException() << "cannot run flat subprocesses with unweighted sampling.";
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
    const StandardEventHandler& eh = *eventHandler();
    const StandardXComb& xc = *eh.xCombs()[b];
    if ( eventHandler()->weighted() && theSamplingBias ) {
      s->enhanceInitialPoints(theSamplingBias->enhanceInitialPoints(xc));
      s->oversamplingFactor(theSamplingBias->oversamplingFactor(xc));
    }
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

  long excptTries = 0;

  gotCrossSections = false;

  lastSampler = samplers.upper_bound(UseRandom::rnd())->second;

  while ( true ) {

    try {
      lastSampler->generate(eventHandler()->weighted());
    } catch (BinSampler::NewMaximum&) {
      continue;
    } catch(BinSampler::UpdateCrossSections) {
      updateCrossSections();
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    } catch (...) {
      throw;
    }

    if ( isnan(lastSampler->lastWeight()) || isinf(lastSampler->lastWeight()) ) {
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    ++theAttempts;

    if ( eventHandler()->weighted() && lastSampler->lastWeight() == 0.0 ) {
      lastSampler->accept();
      lastSampler = samplers.upper_bound(UseRandom::rnd())->second;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    if ( theAlmostUnweighted ) {
      double w = abs(lastSampler->lastWeight()/lastSampler->bias());
      if ( w <= theMaxWeight ) {
	if ( UseRandom::rnd() > w/theMaxWeight ) {
	  if ( ++excptTries == eventHandler()->maxLoop() )
	    break;
	  continue;
	}
      }
    }

    break;

  }

  ++theAccepts;

  if ( excptTries == eventHandler()->maxLoop() )
    throw Exception()
      << "GeneralSampler::generate() : Maximum number of tries to re-run event "
      << "selection reached. Aborting now." << Exception::runerror;

  lastPoint() = lastSampler->lastPoint();
  lastSampler->accept();

  if ( !eventHandler()->weighted() ) {
    theSumWeights += sign(lastSampler->lastWeight());
    theSumWeights2 += 1.0;
    return sign(lastSampler->lastWeight());
  } else {
    double w = lastSampler->lastWeight()/lastSampler->bias();
    if ( theAlmostUnweighted ) {
      if ( w <= theMaxWeight )
	w = theMaxWeight*sign(lastSampler->lastWeight()/lastSampler->bias());
    }
    theSumWeights += w;
    theSumWeights2 += sqr(w);
    return w;
  }
  return 0.;

}

void GeneralSampler::rejectLast() {
  lastSampler->reject();
  if ( !eventHandler()->weighted() ) {
    theSumWeights -= sign(lastSampler->lastWeight());
    theSumWeights2 -= 1.0;
  } else {
    double w = lastSampler->lastWeight()/lastSampler->bias();
    if ( theAlmostUnweighted ) {
      if ( w <= theMaxWeight )
	w = theMaxWeight*sign(lastSampler->lastWeight()/lastSampler->bias());
    }
    theSumWeights -= w;
    theSumWeights2 -= sqr(w);
  }
  --theAttempts;
  --theAccepts;
}

void GeneralSampler::currentCrossSections() const {

  if ( !isSampling )
    return;

  if ( gotCrossSections )
    return;

  if ( crossSectionCalls > 0 ) {
    if ( ++crossSectionCalls == theUpdateAfter ) {
      crossSectionCalls = 0;
    } else return;
  }

  ++crossSectionCalls;

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

  theMaxWeight = 0.0;

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    if ( (isSampling && s->second == lastSampler) ||
	 !isSampling )
      s->second->nextIteration();
    theMaxWeight = max(theMaxWeight,abs(s->second->iterations().back().maxWeight()));
    if ( isSampling && s->second == lastSampler ) {
      s->second->maxWeight(s->second->iterations().back().maxWeight());
      s->second->minWeight(s->second->iterations().back().minWeight());
    }
    xsec += s->second->averageWeight();
    var += s->second->averageWeightVariance();
    if ( !theFlatSubprocesses ) {
      sumbias += s->second->averageAbsWeight() * s->second->oversamplingFactor();
    } else {
      sumbias += 1.;
    }
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

  double minBias = !theFlatSubprocesses ? -1.0 : 1/sumbias;

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers.begin();
	s != samplers.end(); ++s ) {
    double abssw = s->second->averageAbsWeight() * s->second->oversamplingFactor();
    if ( abssw == 0.0 )
      continue;
    if ( !theFlatSubprocesses ) {
      minBias = minBias > 0.0 ? min(minBias,abssw/sumbias) : abssw/sumbias;
      s->second->bias(abssw/sumbias);
      current += abssw;
    } else {
      s->second->bias(1./sumbias);
      current += 1.;
    }
    newSamplers[current/sumbias] = s->second;
  }

  theMaxWeight /= minBias;

  samplers = newSamplers;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void GeneralSampler::dofinish() {
  set<string> compensating;
  for ( map<double,Ptr<BinSampler>::ptr>::const_iterator s =
	  samplers.begin(); s != samplers.end(); ++s ) {
    if ( s->second->compensating() ) {
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

  if ( runCombinationData ) {

    string dataName = generator()->filename() + "-sampling.dat";

    ofstream data(dataName.c_str());

    double runXSec =
      theSumWeights/theAttempts;
    double runXSecErr =
      (1./theAttempts)*(1./(theAttempts-1.))*
      abs(theSumWeights2 - sqr(theSumWeights)/theAttempts);
      
    data << setprecision(20);

    data << "CrossSectionCombined "
	 << (integratedXSec()/nanobarn) << " +/- "
	 << (integratedXSecErr()/nanobarn) << "\n"
	 << "CrossSectionRun "
	 << runXSec << " +/- " << sqrt(runXSecErr) << "\n"
	 << "PointsAttempted " << theAttempts << "\n"
	 << "PointsAccepted " << theAccepts << "\n"
	 << "SumWeights " << theSumWeights << "\n"
	 << "SumWeights2 " << theSumWeights2 << "\n"
	 << flush;

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
  os << theBinSampler << theSamplingBias << theVerbose << theFlatSubprocesses 
     << samplers << lastSampler << theUpdateAfter
     << theIntegratedXSec << theIntegratedXSecErr
     << theSumWeights << theSumWeights2 
     << theAttempts << theAccepts
     << norm << runCombinationData << theMaxWeight
     << theAlmostUnweighted;
}

void GeneralSampler::persistentInput(PersistentIStream & is, int) {
  is >> theBinSampler >> theSamplingBias >> theVerbose >> theFlatSubprocesses 
     >> samplers >> lastSampler >> theUpdateAfter
     >> theIntegratedXSec >> theIntegratedXSecErr
     >> theSumWeights >> theSumWeights2 
     >> theAttempts >> theAccepts
     >> norm >> runCombinationData >> theMaxWeight
     >> theAlmostUnweighted;
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


  static Reference<GeneralSampler,SamplingBias> interfaceSamplingBias
    ("SamplingBias",
     "Set the sampling bias to apply.",
     &GeneralSampler::theSamplingBias, false, false, true, true, false);


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

  static Switch<GeneralSampler,bool> interfaceRunCombinationData
    ("RunCombinationData",
     "",
     &GeneralSampler::runCombinationData, false, false, false);
  static SwitchOption interfaceRunCombinationDataOn
    (interfaceRunCombinationData,
     "On",
     "",
     true);
  static SwitchOption interfaceRunCombinationDataOff
    (interfaceRunCombinationData,
     "Off",
     "",
     false);

  interfaceRunCombinationData.rank(-1);

  static Switch<GeneralSampler,bool> interfaceAlmostUnweighted
    ("AlmostUnweighted",
     "",
     &GeneralSampler::runCombinationData, false, false, false);
  static SwitchOption interfaceAlmostUnweightedOn
    (interfaceAlmostUnweighted,
     "On",
     "",
     true);
  static SwitchOption interfaceAlmostUnweightedOff
    (interfaceAlmostUnweighted,
     "Off",
     "",
     false);

  interfaceAlmostUnweighted.rank(-1);

}

