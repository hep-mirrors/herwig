// -*- C++ -*-
//
// GeneralSampler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

#include "ThePEG/Utilities/ColourOutput.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include "Herwig/API/RunDirectories.h"

#include "Herwig/Utilities/XML/ElementIO.h"

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <sstream> 

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
    maximumExceededBy(0.), correctWeights(0.),theMaxEnhancement(1.05), didReadGrids(false),
    theParallelIntegration(false),
  theIntegratePerJob(0), theIntegrationJobs(0), theIntegrationJobsCreated(0),
  justAfterIntegrate(false), theWriteGridsOnFinish(false) {}

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

  if ( theParallelIntegration &&
       runLevel() == ReadMode )
    throw Exception()
      << "\n----------------------------------------------------\n\n"
      << "Parallel integration is only supported\n in the build/integrate/run mode\n\n"
      << "----------------------------------------------------\n"
      << Exception::abortnow;

  if ( runLevel() == ReadMode ||
       runLevel() == IntegrationMode ) {
    assert(theSamplers.empty());
  }

  if ( theParallelIntegration ) {

    if ( !theIntegratePerJob && !theIntegrationJobs )
      throw Exception()
	<< "Please specify the number of subprocesses per integration job or the "
	<< "number of integration jobs to be created."
	<< Exception::abortnow;

    if ( theIntegrationJobs ) {
      unsigned int nintegrate = eventHandler()->nBins()/theIntegrationJobs;
      if ( eventHandler()->nBins() % theIntegrationJobs != 0 )
	++nintegrate;
      theIntegratePerJob = max(theIntegratePerJob,nintegrate);
    }

    unsigned int jobCount = 0;

    ofstream* jobList = 0;

    generator()->log() 
      << "----------------------------------------------------\n"
      << "preparing integration jobs ...\n" << flush;

    vector<int> randomized;
    vector<int> pickfrom;
    for ( int b = 0; b < eventHandler()->nBins(); ++b )
      pickfrom.push_back(b);

    //set<int> check;

    while ( !pickfrom.empty() ) {
      size_t idx = UseRandom::irnd(pickfrom.size());
      randomized.push_back(pickfrom[idx]);
      pickfrom.erase(pickfrom.begin() + idx);
    }

    int b = 0;
    for ( vector<int>::const_iterator bx = randomized.begin();
	  bx != randomized.end(); ++bx, ++b ) {

      if ( b == 0 || b % theIntegratePerJob == 0 ) {
	if ( jobList ) {
	  jobList->close();
	  delete jobList;
	  jobList = 0;
	}
	ostringstream name;
	string prefix = RunDirectories::runStorage();
	if ( prefix.empty() )
	  prefix = "./";
	else if ( *prefix.rbegin() != '/' )
	  prefix += "/";
	name << prefix << "integrationJob" << jobCount<<"List";
	++jobCount;
	string fname = name.str();
	jobList = new ofstream(fname.c_str());
	if ( !*jobList ) {
	  delete jobList;
	  throw Exception() << "Failed to write integration job list"
			    << Exception::abortnow;
	}
      }

      *jobList << *bx << " ";

    }
    theIntegrationJobsCreated = jobCount;

    generator()->log() 
      << "---------------------------------------------------\n\n"
      << "Wrote " << jobCount << " integration jobs\n"
      << "Please submit integration jobs with the\nintegrate --jobid=x\ncommand for job ids "
      << "from 0 to " << (jobCount-1) << "\n\n"
      << "e.g.:\n\n" << ANSI::yellow
      << " for i in $(seq 0 "<< (jobCount-1) <<");do Herwig integrate --jobid=$i "<<generator()->runName()<<".run & done \n\n" << ANSI::reset
      << "---------------------------------------------------\n"
      << flush;

    if ( jobList ) {
      jobList->close();
      delete jobList;
      jobList = 0;
    }

    theParallelIntegration = false;
    return;

  }

  if ( runLevel() == BuildMode )
    return;

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
  if ( integrationList() != "" ) {
    string prefix = RunDirectories::runStorage();
    assert ( !prefix.empty() );
    string fname = prefix.substr(0, prefix.size()-1) + "List";

    ifstream jobList(fname.c_str());
    if ( jobList ) {
      int b = 0;
      while ( jobList >> b )
	binsToIntegrate.insert(b);
    } else {
      Repository::clog() 
	<< "Job list '"
	<< fname << "' not found.\n"
	<< "Assuming empty integration job\n" << flush;
      return;
    }
  }

  if ( binsToIntegrate.empty() ) {
    for ( int b = 0; b < eventHandler()->nBins(); ++b ) 
      binsToIntegrate.insert(b);
  }

  boost::progress_display* progressBar = 0;
  if ( !theVerbose && !justAfterIntegrate ) {
    Repository::clog() << "integrating subprocesses";
    progressBar = new boost::progress_display(binsToIntegrate.size(),Repository::clog());
  }

  int count=0;
  bool reuseGrid = false;
  bool missingGrid = false;
  for ( set<int>::const_iterator bit = binsToIntegrate.begin(); bit != binsToIntegrate.end(); ++bit ) {
    count++;
    if(theVerbose&&
       (runLevel() == ReadMode ||
        runLevel() == IntegrationMode))
       cout<<"\nIntegrate "<< count <<" of "<<binsToIntegrate.size() <<":\n"<<flush;
    Ptr<BinSampler>::ptr s = theBinSampler->cloneMe();
    s->eventHandler(eventHandler());
    s->sampler(this);
    s->bin(*bit);
    lastSampler(s);
    s->doWeighted(eventHandler()->weighted());
    s->setupRemappers(theVerbose);
    if ( justAfterIntegrate )
      s->readIntegrationData();
    reuseGrid = reuseGrid || s->existsGrid();
    missingGrid = missingGrid || ( ! s->existsGrid() );
    s->initialize(theVerbose);
    samplers()[*bit] = s;
    if ( !theVerbose && !justAfterIntegrate )
      ++(*progressBar);
    if ( s->nanPoints() && theVerbose ) {
      Repository::clog() << "warning: " 
			 << s->nanPoints() << " of "
			 << s->allPoints() << " points with nan or inf weight.\n"
			 << flush;
    }
  }

  if ( progressBar ) {
    delete progressBar;
    progressBar = 0;
  }

  if ( missingGrid && runLevel() == RunMode )
    generator()->log()
    << "\n----------------------------------------------------\n\n"
    << "Warning:No grid file could be found at the start of\n"
    << "this run.\n\n"
    << "* For a read/run setup intented to be used with \n"
    << "  --setupfile \n"
    << "  please consider using the build/integrate/run setup.\n"
    << "* For a build/integrate/run setup to be used with\n"
    << "  --setupfile\n"
    << "  please ensure that the same setupfile is provided\n"
    << "  to both the integrate and run steps.\n\n"
    << "---------------------------------------------------\n" << flush;

  if ( runLevel() == ReadMode ||
       runLevel() == IntegrationMode ) {
    if ( reuseGrid )
      Repository::clog()
        << "----------------------------------------------------\n\n"
        << "Re-using an existing grid as starting point for grid\n"
        << "optimization. Please consider removing the grid files \n"
        << "and re-running the grid adaption when there have\n"
        << "been significant changes to parameters, cuts, etc.\n\n"
        << "----------------------------------------------------\n"
        << flush;
  }

  if ( runLevel() == IntegrationMode ) {
    theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
    for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	  s != samplers().end(); ++s ) {
      s->second->saveGrid();
      s->second->saveRemappers();
      s->second->saveIntegrationData();
    }
    writeGrids();
    return;
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

  if ( !justAfterIntegrate ) {
    theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
    for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	  s != samplers().end(); ++s ) {
      s->second->saveGrid();
      s->second->saveRemappers();
    }
    writeGrids();
  }

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
      else if ( p < UseRandom::rnd() ){
        weight = 0.0;
        // The lastSampler was picked according to the bias of the process.
        --excptTries;
      }
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

    if ( ! isfinite(lastSampler()->lastWeight()) ) {
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
      correctWeights+=weight;
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
    allMax = max(allMax,s->second->maxWeight()*theMaxEnhancement);
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
    double wref = theGlobalMaximumWeight ? allMax : 
                                 s->second->maxWeight()*theMaxEnhancement;
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

void GeneralSampler::prepare() {
  readGrids();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void GeneralSampler::doinit() {
  if ( RunDirectories::empty() )
    RunDirectories::pushRunId(generator()->runName());
  if ( integratePerJob() || integrationJobs() ) {
    theParallelIntegration = true;
    theIntegratePerJob = integratePerJob();
    theIntegrationJobs = integrationJobs();
  }
  readGrids();
  bool missingGrid = false;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) 
    missingGrid = missingGrid || ( ! s->second->existsGrid() );
  if ( missingGrid && runLevel() == RunMode )
    generator()->log()
      << "\n---------------------------------------------------\n\n"
      << "Warning: No grid file could be found at the start of this run.\n\n"
      << "* For a read/run setup intented to be used with --setupfile please consider\n"
      << "  using the build/integrate/run setup.\n"
      << "* For a build/integrate/run setup to be used with --setupfile please ensure\n"
      << "  that the same setupfile is provided to both, the integrate and run steps.\n\n"
      << "---------------------------------------------------\n" << flush;
  if ( samplers().empty() && runLevel() == RunMode )
    justAfterIntegrate = true;
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
    //generator()->log() << maximumExceeds << " of " << theAttempts
		//	    << " attempted points exceeded the guessed maximum weight\n"
		//	    << "with an average relative deviation of "
		//       << maximumExceededBy/maximumExceeds << "\n\n" << flush;
    generator()->log() <<"\n\n\nNote: In this run "<<maximumExceeds<<" of the "<<theAccepts<<" accepted events\n"
                       <<"were found with a weight W larger than the expected Wmax.\n";
    
generator()->log() <<"This corresponds to a cross section difference between:\n"
		<<"   UnitWeights:       "<< theMaxWeight*theSumWeights/theAttempts<<"nb\n"
                       <<"   AlmostUnweighted:  "<< theMaxWeight*correctWeights/theAttempts<< "nb\n"<<
                      " use 'set Sampler:AlmostUnweighted Yes' to switch to non-unit weights.\n\n";

    generator()->log() <<"The maximum weight determined in the read/integrate step has been enhanced by \n"<<
                         "   set /Herwig/Samplers/Sampler:MaxEnhancement "<< theMaxEnhancement<<
                         ".\nIf the rate of excessions ("<<(double)maximumExceeds*100/(double)theAccepts<<
                         "%) or the change of the cross section is large,\nyou can try to:\n\n"<<
                         "Enhance the number of points used in the read/integrate step\n"<<
                         "   set /Herwig/Samplers/Sampler:BinSampler:InitialPoints ...\n\n"<<
                         "and/or enhance the reference weight found in the read/integrate step\n"<<  
                         "   set /Herwig/Samplers/Sampler:MaxEnhancement 1.x\n\n"<<
                         "If this does not help (and your process is well defined by cuts)\n"<<
                         "don't hesitate to contact herwig@projects.hepforge.org.\n\n";

  }

  if ( runCombinationData ) {

    string dataName = RunDirectories::runStorage();
    if ( dataName.empty() )
      dataName = "./";
    else if ( *dataName.rbegin() != '/' )
      dataName += "/";
    dataName += "HerwigSampling.dat";

    ofstream data(dataName.c_str());

    double runXSec =
      theMaxWeight*theSumWeights/theAttempts;
    double runXSecErr =
      sqr(theMaxWeight)*(1./theAttempts)*(1./(theAttempts-1.))*
      abs(theSumWeights2 - sqr(theSumWeights)/theAttempts);
      
    data << setprecision(17);

    data << "CrossSectionCombined "
	 << (integratedXSec()/nanobarn) << " +/- "
	 << (integratedXSecErr()/nanobarn) << "\n"
	 << "CrossSectionRun "
	 << runXSec << " +/- " << sqrt(runXSecErr) << "\n"
	 << "PointsAttempted " << theAttempts << "\n"
	 << "PointsAccepted " << theAccepts << "\n"
	 << "SumWeights " << theSumWeights*theMaxWeight << "\n"
	 << "SumWeights2 " << theSumWeights2*sqr(theMaxWeight) << "\n"
	 << flush;

  }

  theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) {
    s->second->saveGrid();
    s->second->saveRemappers();
    if ( justAfterIntegrate )
      s->second->saveIntegrationData();
  }

  if ( theWriteGridsOnFinish )
    writeGrids();

  SamplerBase::dofinish();

}

void GeneralSampler::doinitrun() {
  readGrids();

  if ( samplers().empty() ) {
    justAfterIntegrate = true;
    if ( !hasSetupFile() )
      initialize();
  } else {
    for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	  s != samplers().end(); ++s ) {
      s->second->setupRemappers(theVerbose);
      if ( justAfterIntegrate )
	s->second->readIntegrationData();
      s->second->initialize(theVerbose);
    }
  }

  bool missingGrid = false;
  for ( map<double,Ptr<BinSampler>::ptr>::iterator s = samplers().begin();
	s != samplers().end(); ++s ) 
    missingGrid = missingGrid || ( ! s->second->existsGrid() );
  if ( missingGrid && !didReadGrids )
    generator()->log()
      << "\n----------------------------------------------------\n\n"
      << "Warning:No grid file could be found at the start of\n"
      << "this run.\n\n"
      << "* For a read/run setup intented to be used with \n"
      << "  --setupfile \n"
      << "  please consider using the build/integrate/run setup.\n"
      << "* For a build/integrate/run setup to be used with\n"
      << "  --setupfile\n"
      << "  please ensure that the same setupfile is provided\n"
      << "  to both the integrate and run steps.\n\n"
      << "---------------------------------------------------\n" << flush;

  isSampling = true;
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
  string dataName = RunDirectories::runStorage();
  if ( dataName.empty() )
    dataName = "./";
  else if ( *dataName.rbegin() != '/' )
    dataName += "/";
  dataName += "HerwigGrids.xml";
  ofstream out(dataName.c_str());
  XML::ElementIO::put(theGrids,out);
}

void GeneralSampler::readGrids() {
  // return if grids were already read
  if ( didReadGrids )
    return;
  
  // check for global HerwigGrids.xml file or combine integration jobs to a global HerwigGrids.xml file
  // Show messages of integration job combination only in the first run  (if no global HerwigGrids.xml file is found in one of the directories)
  // or in case of an error
  // Check if a globalHerwigGridsFileFound was found and keep messages in a stringstream buffer beforehand
  bool globalHerwigGridsFileFound = false;
  bool integrationJobCombinationSuccessful = true;
  std::stringstream messageBuffer;  
  
  RunDirectories directories;
  while ( directories && !didReadGrids ) {
    string dataName = directories.nextRunStorage();
    if ( dataName.empty() )
      dataName = "./";
    else if ( *dataName.rbegin() != '/' )
      dataName += "/";
    string directoryName = dataName;
    dataName += "HerwigGrids.xml";
    ifstream in(dataName.c_str());
    if ( in ) {
      theGrids = XML::ElementIO::get(in);
      didReadGrids = true;
      // Set to true if in any of the directories a global HerwigGrid.xml file was found
      globalHerwigGridsFileFound = true;
    }
    else {
      // Check if integrationJob was split and try to merge single integrationJobs together
      // integrationJobsCreated() == 0 indicates that parallel integration has not been
      // requested, while the parallel integration parameters may well yield a single job
      if(integrationJobsCreated() >= 1 && runLevel() == RunMode) {
	messageBuffer << "\n\n* Global HerwigGrids.xml file does not exist yet"
		      << "\n  and integration jobs were split into " << integrationJobsCreated() << " integration jobs."
		      << "\n  Trying to combine single integration jobs to a global HerwigGrids.xml file"
		      << "\n  using the following directory " << directoryName  << ".";


	theGrids = XML::Element(XML::ElementTypes::Element,"Grids");

	integrationJobCombinationSuccessful = true;
			    
	for(unsigned int currentProcessedIntegrationJobNum = 0; currentProcessedIntegrationJobNum < integrationJobsCreated(); ++currentProcessedIntegrationJobNum) {
    ostringstream currentProcessedIntegrationJob;
	  currentProcessedIntegrationJob << directoryName << "integrationJob" << currentProcessedIntegrationJobNum << "/HerwigGrids.xml";
	  if(boost::filesystem::exists(boost::filesystem::path(currentProcessedIntegrationJob.str()))) {
	    ifstream localGridFileIN(currentProcessedIntegrationJob.str().c_str());
	    if(localGridFileIN) {
	      theGrids = theGrids + XML::ElementIO::get(localGridFileIN);
	      messageBuffer  << "\n* Added integration job " << currentProcessedIntegrationJobNum << " to global HerwigGrids.xml file.";
	    } 
	    else {
	      integrationJobCombinationSuccessful = false;
	      messageBuffer << "\n* Could not open/add integration job " << currentProcessedIntegrationJobNum << " to global HerwigGrids.xml file.";
	    }
	  }  
	  else {
	    integrationJobCombinationSuccessful = false;
	    messageBuffer << "\n* Could not find integration job " << currentProcessedIntegrationJob.str();
	  }
	}
	
	if(integrationJobCombinationSuccessful) {
	  string globalGridFile = directoryName + "HerwigGrids.xml";
	  ofstream globalGridFileOF(globalGridFile.c_str());
	  XML::ElementIO::put(theGrids,globalGridFileOF);
	  messageBuffer << "\n* Global HerwigGrids.xml file was created, the integration jobs 0 to " << integrationJobsCreated()-1 
			<< " were combined."
			<< "\n* If previous warnings in regards to the HerwigGrids.xml file occured, these can be safely ignored."
			<< "\n* Note: This message will occur only in the first run and will be suppressed in further runs.\n" 
			<< flush;
	  didReadGrids = true;
	}
	else {
	  messageBuffer << "\n* Global HerwigGrids.xml file could not be created due to failed combination of integration jobs."
				  << "\n  Please check the above-mentioned missing/failed integration jobs which are needed for the combination."
				  << "\n* Note: It can be that the HerwigGrids.xml file is searched and can be found in further directories."
				  << "\n  In this case you can ignore this warning message.\n" << flush;
	}
      }  
    }  
  }
  
  // Show messages if global HerwigGrids.xml file was not found or first combination run 
  if (!globalHerwigGridsFileFound && (theVerbose || !integrationJobCombinationSuccessful))
    BaseRepository::cout() << messageBuffer.str() << "\n" << flush;
  
  if ( !didReadGrids )
    theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
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
     << maximumExceededBy << correctWeights << theMaxEnhancement     
     << theParallelIntegration
     << theIntegratePerJob << theIntegrationJobs 
     << theIntegrationJobsCreated << theWriteGridsOnFinish;
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
     >> maximumExceededBy >> correctWeights >> theMaxEnhancement
     >> theParallelIntegration
     >> theIntegratePerJob >> theIntegrationJobs 
     >> theIntegrationJobsCreated >> theWriteGridsOnFinish;
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
  static SwitchOption interfaceVerboseYes
    (interfaceVerbose,
     "Yes",
     "",
     true);
  static SwitchOption interfaceVerboseNo
    (interfaceVerbose,
     "No",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceAddUpSamplers
    ("AddUpSamplers",
     "Calculate cross sections from adding up individual samplers.",
     &GeneralSampler::theAddUpSamplers, false, false, false);
  static SwitchOption interfaceAddUpSamplersYes
    (interfaceAddUpSamplers,
     "Yes",
     "",
     true);
  static SwitchOption interfaceAddUpSamplersNo
    (interfaceAddUpSamplers,
     "No",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceGlobalMaximumWeight
    ("GlobalMaximumWeight",
     "Use a global maximum weight instead of partial unweighting.",
     &GeneralSampler::theGlobalMaximumWeight, true, false, false);
  static SwitchOption interfaceGlobalMaximumWeightYes
    (interfaceGlobalMaximumWeight,
     "Yes",
     "",
     true);
  static SwitchOption interfaceGlobalMaximumWeightNo
    (interfaceGlobalMaximumWeight,
     "No",
     "",
     false);


  static Parameter<GeneralSampler,double> interfaceMaxEnhancement
    ("MaxEnhancement",
     "Enhance the maximum reference weight found in the read step.",
     &GeneralSampler::theMaxEnhancement, 1.1, 1.0, 1.5,
     false, false, Interface::limited);


  static Switch<GeneralSampler,bool> interfaceFlatSubprocesses
    ("FlatSubprocesses",
     "[debug] Perform a flat subprocess selection.",
     &GeneralSampler::theFlatSubprocesses, false, false, false);
  static SwitchOption interfaceFlatSubprocessesYes
    (interfaceFlatSubprocesses,
     "Yes",
     "",
     true);
  static SwitchOption interfaceFlatSubprocessesNo
    (interfaceFlatSubprocesses,
     "No",
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
  static SwitchOption interfaceRunCombinationDataYes
    (interfaceRunCombinationData,
     "Yes",
     "",
     true);
  static SwitchOption interfaceRunCombinationDataNo
    (interfaceRunCombinationData,
     "No",
     "",
     false);

  static Switch<GeneralSampler,bool> interfaceAlmostUnweighted
    ("AlmostUnweighted",
     "",
     &GeneralSampler::theAlmostUnweighted, false, false, false);
  static SwitchOption interfaceAlmostUnweightedYes
    (interfaceAlmostUnweighted,
     "Yes",
     "",
     true);
  static SwitchOption interfaceAlmostUnweightedNo
    (interfaceAlmostUnweighted,
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

  static Parameter<GeneralSampler,unsigned int> interfaceIntegrationJobs
    ("IntegrationJobs",
     "The maximum number of integration jobs to create.",
     &GeneralSampler::theIntegrationJobs, 0, 0, 0,
     false, false, Interface::lowerlim);
    
  static Parameter<GeneralSampler,unsigned int> interfaceIntegrationJobsCreated
    ("IntegrationJobsCreated",
     "The number of integration jobs which were actually created.",
     &GeneralSampler::theIntegrationJobsCreated, 1, 1, 0,
     false, false, Interface::lowerlim);    

  static Switch<GeneralSampler,bool> interfaceWriteGridsOnFinish
    ("WriteGridsOnFinish",
     "Write grids on finishing a run.",
     &GeneralSampler::theWriteGridsOnFinish, false, false, false);
  static SwitchOption interfaceWriteGridsOnFinishYes
    (interfaceWriteGridsOnFinish,
     "Yes",
     "",
     true);
  static SwitchOption interfaceWriteGridsOnFinishNo
    (interfaceWriteGridsOnFinish,
     "No",
     "",
     false);

}

