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
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/LoopGuard.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include "Herwig++/Utilities/XML/ElementIO.h"

#include <boost/progress.hpp>
#include <cstdlib>

using namespace Herwig;

GeneralSampler::GeneralSampler() 
  : theVerbose(false),
    theIntegratedXSec(ZERO), theIntegratedXSecErr(ZERO),
    theUpdateAfter(1), crossSectionCalls(0), gotCrossSections(false),
    theSumWeights(0.), theSumWeights2(0.), 
    theAttempts(0), theAccepts(0),
    theMaxWeight(0.0), theAddUpSamplers(false),
    theGlobalMaximumWeight(true), theFlatSubprocesses(false),
    isSampling(false), theMinSelection(0.01), runCombinationData(false),
    theAlmostUnweighted(false), maximumExceeds(0),
    maximumExceededBy(0.), didReadGrids(false), thePostponeInitialize(false),
    theParallelIntegration(false), theSaveStatistics(false),
    theIntegratePerJob(0), theIgnoreIntegrationData(true) {}

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

  if ( theParallelIntegration ) {

    if ( !theIntegratePerJob )
      throw Exception()
	<< "Please specify the number of integration jobs"
	<< Exception::abortnow;

    unsigned int jobCount = 0;

    ofstream* jobList = 0;

    for ( int b = 0; b < eventHandler()->nBins(); ++b ) {

      if ( b == 0 || (b+1) % theIntegratePerJob == 0 ) {
	if ( jobList ) {
	  jobList->close();
	  jobList = 0;
	}
	ostringstream name;
	name << "./HerwigIntegration/run" << jobCount;
	++jobCount;
	string cmd = "mkdir -p " + name.str();
	std::system(cmd.c_str());
	name << "/integrationBins.dat";
	string fname = name.str();
	jobList = new ofstream(fname.c_str());
	if ( !*jobList )
	  throw Exception() << "Failed to write integratino job list"
			    << Exception::abortnow;
      }

      *jobList << b << " ";

    }

    if ( jobList ) {
      jobList->close();
      jobList = 0;
    }

    theParallelIntegration = false;
    theSaveStatistics = true;
    if ( postponeInitialize() )
      thePostponeInitialize = false;
    return;

  }

  if ( postponeInitialize() ) {
    thePostponeInitialize = false;
    return;
  }

  if ( !samplers().empty() )
    return;

  if ( binSampler()->adaptsOnTheFly() ) {
    if ( !theAddUpSamplers ) {
      Repository::clog() << "Warning: On-the-fly adapting samplers require cross section calculation from "
			 << "adding up individual samplers. The AddUpSamplers flag will be switched on.";
    }
    theAddUpSamplers = true;
  }

  if ( !weighted() && !binSampler()->canUnweight() )
    throw Exception() << "Unweighted events requested from weighted bin sampler object.";

  if ( theFlatSubprocesses && !theGlobalMaximumWeight ) {
      Repository::clog() << "Warning: Can only use a global maximum weight when selecting subprocesses "
			 << "uniformly. The GlobalMaximumWeight flag will be switched on.";
    theGlobalMaximumWeight = true;
  }

  set<int> binsToIntegrate;
  bool integrationJob = false;
  if ( theIntegratePerJob ) {
    ifstream jobList("integrationBins.dat");
    if ( !jobList )
      throw Exception() << "Failed to load integration job list"
			<< Exception::abortnow;
    int b = 0;
    while ( jobList >> b )
      binsToIntegrate.insert(b);
    integrationJob = true;
  }

  if ( binsToIntegrate.empty() ) {
    for ( int b = 0; b < eventHandler()->nBins(); ++b ) 
      binsToIntegrate.insert(b);
  }

  boost::progress_display* progressBar = 0;
  if ( !theVerbose ) {
    Repository::clog() << "integrating subprocesses";
    progressBar = new boost::progress_display(binsToIntegrate.size(),Repository::clog());
  }

  for ( set<int>::const_iterator bit = binsToIntegrate.begin(); bit != binsToIntegrate.end(); ++bit ) {
    Ptr<BinSampler>::ptr s = theBinSampler->cloneMe();
    s->eventHandler(eventHandler());
    s->sampler(this);
    s->bin(*bit);
    lastSampler(s);
    s->doWeighted(eventHandler()->weighted());
    s->setupRemappers(theVerbose);
    if ( !theIgnoreIntegrationData &&
	 !integrationJob )
      s->readIntegrationData();
    s->initialize(theVerbose);
    samplers()[*bit] = s;
    if ( !theVerbose )
      ++(*progressBar);
    if ( s->nanPoints() && theVerbose ) {
      Repository::clog() << "warning: " 
			 << s->nanPoints() << " of "
			 << s->allPoints() << " points with nan or inf weight.\n"
			 << flush;
    }
  }

  if ( theVerbose ) {
    bool oldAdd = theAddUpSamplers;
    theAddUpSamplers = true;
    try {
      Repository::clog() << "estimated total cross section is ( "
			 << integratedXSec()/nanobarn << " +/- "
			 << integratedXSecErr()/nanobarn << " ) nb\n" << flush;
    } catch (...) {
      theAddUpSamplers = oldAdd;
      throw;
    }
    theAddUpSamplers = oldAdd;
  }

  updateSamplers();

  if ( samplers().empty() ) {
    throw Exception() << "No processes with non-zero cross section present."
		      << Exception::abortnow;
  }

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    s->second->saveRemappers();
    if ( theSaveStatistics )
      s->second->saveIntegrationData();
  }

  writeGrids();

}

double GeneralSampler::generate() {

  long excptTries = 0;

  gotCrossSections = false;

  lastSampler(samplers().upper_bound(UseRandom::rnd())->second);

  double weight = 0.;

  while ( true ) {

    try {
      weight = 1.0;
      double p = lastSampler()->referenceWeight()/lastSampler()->bias()/theMaxWeight;
      if ( weighted() )
	weight *= p;
      else if ( p < UseRandom::rnd() )
	weight = 0.0;
      if ( weight != 0.0 )
	weight *= lastSampler()->generate()/lastSampler()->referenceWeight();
    } catch(BinSampler::NextIteration) {
      updateSamplers();
      lastSampler(samplers().upper_bound(UseRandom::rnd())->second);
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    } catch (...) {
      throw;
    }

    if ( isnan(lastSampler()->lastWeight()) || isinf(lastSampler()->lastWeight()) ) {
      lastSampler() = samplers().upper_bound(UseRandom::rnd())->second;
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    theAttempts += 1;

    if ( abs(weight) == 0.0 ) {
      lastSampler(samplers().upper_bound(UseRandom::rnd())->second);
      if ( ++excptTries == eventHandler()->maxLoop() )
	break;
      continue;
    }

    if ( !eventHandler()->weighted() && !theAlmostUnweighted ) {
      if ( abs(weight) > 1. ) {
	++maximumExceeds;
	maximumExceededBy += abs(weight)-1.;
      }
      if ( weight > 0.0 )
	weight = 1.;
      else
	weight = -1.;
    }

    break;

  }

  theAccepts += 1;

  if ( excptTries == eventHandler()->maxLoop() )
    throw Exception()
      << "GeneralSampler::generate() : Maximum number of tries to re-run event "
      << "selection reached. Aborting now." << Exception::runerror;

  lastPoint() = lastSampler()->lastPoint();
  lastSampler()->accept();

  theSumWeights += weight;
  theSumWeights2 += sqr(weight);

  return weight;

}

void GeneralSampler::rejectLast() {
  if ( !lastSampler() )
    return;
  double w = 0.0;
  if ( weighted() )
    w = lastSampler()->lastWeight()/lastSampler()->bias()/theMaxWeight;
  else
    w = lastSampler()->lastWeight()/lastSampler()->referenceWeight();
  lastSampler()->reject();
  theSumWeights -= w;
  theSumWeights2 -= sqr(w);
  theAttempts -= 1;
  theAccepts -= 1;
}

void GeneralSampler::updateSamplers() {

  map<double,Ptr<BinSampler>::ptr> checkedSamplers;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    if ( s->second->averageAbsWeight() == 0.0 ) {
      generator()->log() << "Warning: no phase space points with non-zero cross section\n"
			 << "could be obtained for the process: "
			 << s->second->process() << "\n"
			 << "This process will not be considered. Try increasing InitialPoints.\n"
			 << flush;
      if ( s->second->nanPoints() ) {
	generator()->log() << "Warning: " 
			   << s->second->nanPoints() << " of "
			   << s->second->allPoints() << " points with nan or inf weight\n"
			   << "in " << s->second->process() << "\n" << flush;
      }
      continue;
    }
    checkedSamplers.insert(*s);
  }

  theSamplers = checkedSamplers;

  if ( samplers().empty() )
    return;

  double allMax = 0.0;

  double sumbias = 0.;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    double bias = 1.;
    if ( !theFlatSubprocesses )
      bias *= s->second->averageAbsWeight();
    s->second->bias(bias);
    sumbias += bias;
    allMax = max(allMax,s->second->maxWeight());
  }

  double nsumbias = 0.0;
  bool needAdjust = false;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    needAdjust |= s->second->bias()/sumbias < theMinSelection;
    s->second->bias(max(s->second->bias()/sumbias,theMinSelection));
    nsumbias += s->second->bias();
  }

  if ( nsumbias == 0.0 ) {
    samplers().clear();
    return;
  }

  if ( needAdjust ) {
    for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	  s != samplers().end(); ++s ) {
      s->second->bias(s->second->bias()/nsumbias);
    }
  }

  theMaxWeight = 0.0;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    double wref = theGlobalMaximumWeight ? allMax : s->second->maxWeight();
    s->second->referenceWeight(wref);
    theMaxWeight = max(theMaxWeight,wref/s->second->bias());
    if ( (isSampling && s->second == lastSampler()) ||
	 !isSampling )
      s->second->nextIteration();
  }

  map<double,Ptr<BinSampler>::ptr> newSamplers;
  double current = 0.;

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    if ( s->second->bias() == 0.0 )
      continue;
    current += s->second->bias();
    newSamplers[current] = s->second;
  }

  samplers() = newSamplers;

}

void GeneralSampler::currentCrossSections() const {

  if ( !theAddUpSamplers ) {
    double n = attempts();
    if ( n > 1 ) {
      theIntegratedXSec = sumWeights()*maxXSec()/attempts();
      double sw = sumWeights(); double sw2 = sumWeights2();
      theIntegratedXSecErr = maxXSec()*sqrt(abs(sw2/n-sqr(sw/n))/(n-1));
    } else {
      theIntegratedXSec = ZERO;
      theIntegratedXSecErr = ZERO;
    }
    return;
  }

  if ( gotCrossSections )
    return;

  if ( crossSectionCalls > 0 ) {
    if ( ++crossSectionCalls == theUpdateAfter ) {
      crossSectionCalls = 0;
    } else return;
  }

  ++crossSectionCalls;
  gotCrossSections = true;

  theIntegratedXSec = ZERO;
  double var = 0.0;

  for ( map<double,Ptr<BinSampler>::ptr>::const_iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    theIntegratedXSec += s->second->integratedXSec();
    var += sqr(s->second->integratedXSecErr()/nanobarn);
  }

  theIntegratedXSecErr = sqrt(var)*nanobarn;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void GeneralSampler::doinit() {
  readGrids();
  SamplerBase::doinit();
}

void GeneralSampler::dofinish() {

  set<string> compensating;
  for ( map<double,Ptr<BinSampler>::ptr>::const_iterator s =
	  samplers().begin(); s != samplers().end(); ++s ) {
    if ( s->second->compensating() ) {
      compensating.insert(s->second->process());
    }
    if ( s->second->nanPoints() ) {
      generator()->log() << "warning: " 
			 << s->second->nanPoints() << " of "
			 << s->second->allPoints() << " points with nan or inf weight\n"
			 << "in " << s->second->process() << "\n" << flush;
    }
    s->second->finalize(theVerbose);
  }
  if ( theVerbose ) {
    if ( !compensating.empty() ) {
      generator()->log() << "warning: sampling for the following processes is still compensating:\n";
      for ( set<string>::const_iterator c = compensating.begin();
	    c != compensating.end(); ++c )
	generator()->log() << *c << "\n";
    }
    generator()->log() << "final integrated cross section is ( "
		       << integratedXSec()/nanobarn << " +/- "
		       << integratedXSecErr()/nanobarn << " ) nb\n" << flush;
  }
  if ( !compensating.empty() ) {
    generator()->log() << "Warning: Some samplers are still in compensating mode.\n" << flush;
  }
  if ( maximumExceeds != 0 ) {
    generator()->log() << maximumExceeds << " of " << theAttempts
			    << " attempted points exceeded the guessed maximum weight\n"
			    << "with an average relative deviation of "
		       << maximumExceededBy/maximumExceeds << "\n" << flush;
  }

  if ( runCombinationData ) {

    string dataName = generator()->filename() + "-sampling.dat";

    ofstream data(dataName.c_str());

    double runXSec =
      theMaxWeight*theSumWeights/theAttempts;
    double runXSecErr =
      sqr(theMaxWeight)*(1./theAttempts)*(1./(theAttempts-1.))*
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

  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    s->second->saveRemappers();
    if ( theSaveStatistics )
      s->second->saveIntegrationData();
  }

  writeGrids();

  SamplerBase::dofinish();

}

void GeneralSampler::doinitrun() {
  readGrids();
  isSampling = true;
  eventHandler()->initrun();
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    s->second->setupRemappers(theVerbose);
    if ( !theIgnoreIntegrationData )
      s->second->readIntegrationData();
    s->second->initialize(theVerbose);
  }
  SamplerBase::doinitrun();
}

void GeneralSampler::rebind(const TranslationMap & trans) {
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = 
	  samplers().begin(); s != samplers().end(); ++s )
    s->second = trans.translate(s->second);
  SamplerBase::rebind(trans);
}

IVector GeneralSampler::getReferences() {
  IVector ret = SamplerBase::getReferences();
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = 
	  samplers().begin(); s != samplers().end(); ++s )
    ret.push_back(s->second);
  return ret;
}

void GeneralSampler::writeGrids() const {
  if ( theGrids.children().empty() )
    return;
  string dataName = generator()->filename() + "-grids.xml";
  ofstream out(dataName.c_str());
  XML::ElementIO::put(theGrids,out);
}

void GeneralSampler::readGrids() {
  if ( didReadGrids )
    return;
  string dataName = generator()->filename() + "-grids.xml";
  ifstream in(dataName.c_str());
  if ( !in ) {
    theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
    return;
  }
  theGrids = XML::ElementIO::get(in);
  didReadGrids = true;
}

void GeneralSampler::persistentOutput(PersistentOStream & os) const {
  os << theVerbose << theBinSampler << theSamplers << theLastSampler
     << theUpdateAfter << crossSectionCalls << gotCrossSections
     << ounit(theIntegratedXSec,nanobarn) 
     << ounit(theIntegratedXSecErr,nanobarn)
     << theSumWeights << theSumWeights2
     << theAttempts << theAccepts << theMaxWeight
     << theAddUpSamplers << theGlobalMaximumWeight
     << theFlatSubprocesses << isSampling << theMinSelection
     << runCombinationData << theAlmostUnweighted << maximumExceeds
     << maximumExceededBy << thePostponeInitialize
     << theParallelIntegration << theSaveStatistics
     << theIntegratePerJob << theIgnoreIntegrationData;
}

void GeneralSampler::persistentInput(PersistentIStream & is, int) {
  is >> theVerbose >> theBinSampler >> theSamplers >> theLastSampler
     >> theUpdateAfter >> crossSectionCalls >> gotCrossSections
     >> iunit(theIntegratedXSec,nanobarn) 
     >> iunit(theIntegratedXSecErr,nanobarn)
     >> theSumWeights >> theSumWeights2
     >> theAttempts >> theAccepts >> theMaxWeight
     >> theAddUpSamplers >> theGlobalMaximumWeight
     >> theFlatSubprocesses >> isSampling >> theMinSelection
     >> runCombinationData >> theAlmostUnweighted >> maximumExceeds
     >> maximumExceededBy >> thePostponeInitialize
     >> theParallelIntegration >> theSaveStatistics
     >> theIntegratePerJob >> theIgnoreIntegrationData;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<GeneralSampler,SamplerBase>
  describeHerwigGeneralSampler("Herwig::GeneralSampler", "HwSampling.so");

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

  static Switch<GeneralSampler,bool> interfaceAddUpSamplers
    ("AddUpSamplers",
     "Calculate cross sections from adding up individual samplers.",
     &GeneralSampler::theAddUpSamplers, false, false, false);
  static SwitchOption interfaceAddUpSamplersOn
    (interfaceAddUpSamplers,
     "On",
     "",
     true);
  static SwitchOption interfaceAddUpSamplersOff
    (interfaceAddUpSamplers,
     "Off",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceGlobalMaximumWeight
    ("GlobalMaximumWeight",
     "Use a global maximum weight instead of partial unweighting.",
     &GeneralSampler::theGlobalMaximumWeight, true, false, false);
  static SwitchOption interfaceGlobalMaximumWeightOn
    (interfaceGlobalMaximumWeight,
     "On",
     "",
     true);
  static SwitchOption interfaceGlobalMaximumWeightOff
    (interfaceGlobalMaximumWeight,
     "Off",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceFlatSubprocesses
    ("FlatSubprocesses",
     "[debug] Perform a flat subprocess selection.",
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

  static Parameter<GeneralSampler,double> interfaceMinSelection
    ("MinSelection",
     "A minimum subprocess selection probability.",
     &GeneralSampler::theMinSelection, 0.01, 0.0, 1.0,
     false, false, Interface::limited);

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

  static Switch<GeneralSampler,bool> interfaceAlmostUnweighted
    ("AlmostUnweighted",
     "",
     &GeneralSampler::theAlmostUnweighted, false, false, false);
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

  static Switch<GeneralSampler,bool> interfacePostponeInitialize
    ("PostponeInitialize",
     "Postpone initialization to happen after event generator modifications.",
     &GeneralSampler::thePostponeInitialize, false, false, false);
  static SwitchOption interfacePostponeInitializeYes
    (interfacePostponeInitialize,
     "Yes",
     "",
     true);
  static SwitchOption interfacePostponeInitializeNo
    (interfacePostponeInitialize,
     "No",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceParallelIntegration
    ("ParallelIntegration",
     "Prepare parallel jobs for integration.",
     &GeneralSampler::theParallelIntegration, false, false, false);
  static SwitchOption interfaceParallelIntegrationYes
    (interfaceParallelIntegration,
     "Yes",
     "",
     true);
  static SwitchOption interfaceParallelIntegrationNo
    (interfaceParallelIntegration,
     "No",
     "",
     false);

  static Parameter<GeneralSampler,unsigned int> interfaceIntegratePerJob
    ("IntegratePerJob",
     "The number of subprocesses to integrate per job.",
     &GeneralSampler::theIntegratePerJob, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<GeneralSampler,bool> interfaceSaveStatistics
    ("SaveStatistics",
     "",
     &GeneralSampler::theSaveStatistics, false, false, false);
  static SwitchOption interfaceSaveStatisticsYes
    (interfaceSaveStatistics,
     "Yes",
     "",
     true);
  static SwitchOption interfaceSaveStatisticsNo
    (interfaceSaveStatistics,
     "No",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceIgnoreIntegrationData
    ("IgnoreIntegrationData",
     "",
     &GeneralSampler::theIgnoreIntegrationData, true, false, false);
  static SwitchOption interfaceIgnoreIntegrationDataYes
    (interfaceIgnoreIntegrationData,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIgnoreIntegrationDataNo
    (interfaceIgnoreIntegrationData,
     "No",
     "",
     false);

}

