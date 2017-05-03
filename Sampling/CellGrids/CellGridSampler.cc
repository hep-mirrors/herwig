// -*- C++ -*-
//
// CellGridSampler.cpp is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CellGridSampler class.
//

#include "CellGridSampler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <boost/progress.hpp>

#include "CellGridSampler.h"
#include "Herwig/Sampling/GeneralSampler.h"

using namespace Herwig;
using namespace ExSample;

CellGridSampler::CellGridSampler() 
  : BinSampler(), SimpleCellGrid(),
    theExplorationPoints(1000), theExplorationSteps(8),
    theGain(0.3), theEpsilon(0.01),
    theMinimumSelection(0.0001), theLuminositySplits(0),
    theChannelSplits(0), theAllChannelSplits(false),
    theUnweightCells(true) {}

CellGridSampler::~CellGridSampler() {}

IBPtr CellGridSampler::clone() const {
  return new_ptr(*this);
}

IBPtr CellGridSampler::fullclone() const {
  return new_ptr(*this);
}

double CellGridSampler::generate() {
  UseRandom rnd;
  double w = SimpleCellGrid::sample(rnd,*this,lastPoint(),
				    !weighted() && initialized() && theUnweightCells, 
				    !initialized());
  if ( !weighted() && initialized() ) {
    double p = min(abs(w),kappa()*referenceWeight())/(kappa()*referenceWeight());
    double sign = w >= 0. ? 1. : -1.;
    if ( p < 1 && UseRandom::rnd() > p )
      w = 0.;
    else
      w = sign*max(abs(w),referenceWeight()*kappa());
  }
  select(w);
  if ( w != 0.0 )
    accept();
  assert(kappa()==1.||sampler()->almostUnweighted());
  return w;
}

void CellGridSampler::adapt() {
  UseRandom rnd;
  set<SimpleCellGrid*> newCells;
  SimpleCellGrid::adapt(theGain,theEpsilon,newCells);
  SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells,Repository::clog());
  SimpleCellGrid::setWeights();
  SimpleCellGrid::updateIntegral();
  SimpleCellGrid::minimumSelection(theMinimumSelection);
}

void CellGridSampler::saveGrid() const {
  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);
}

bool CellGridSampler::existsGrid() const {
  list<XML::Element>::iterator git = sampler()->grids().children().begin();
  for ( ; git != sampler()->grids().children().end(); ++git ) {
    if ( git->type() != XML::ElementTypes::Element )
      continue;
    if ( git->name() != "CellGrid" )
      continue;
    string proc;
    git->getFromAttribute("process",proc);
    if ( proc == id() ) 
      return true;
  }
  return false;
}

void CellGridSampler::initialize(bool progress) {

  bool haveGrid = false;
  list<XML::Element>::iterator git = sampler()->grids().children().begin();
  for ( ; git != sampler()->grids().children().end(); ++git ) {
    if ( git->type() != XML::ElementTypes::Element )
      continue;
    if ( git->name() != "CellGrid" )
      continue;
    string proc;
    git->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveGrid = true;
      break;
    }
  }

  if ( haveGrid ) {
    SimpleCellGrid::fromXML(*git);
    sampler()->grids().erase(git);
    didReadGrids();
  }

  lastPoint().resize(dimension());
  if (randomNumberString()!="") 
  for(size_t i=0;i<lastPoint().size();i++){
     RandomNumberHistograms[RandomNumberIndex(id(),i)] = make_pair( RandomNumberHistogram(),0.);
  }

  if ( initialized() ) {
    if ( !hasGrids() )
      throw Exception() << "CellGridSampler: Require existing grid when starting to run.\n"
			<< "Did you miss setting --setupfile?"
			<< Exception::abortnow;
    return;
  }

  if ( haveGrid ) {
    if ( !integrated() )
      runIteration(initialPoints(),progress);
    isInitialized();
    return;
  }

  SimpleCellGrid::boundaries(vector<double>(dimension(),0.0),vector<double>(dimension(),1.0));
  SimpleCellGrid::weightInformation().resize(dimension());

  UseRandom rnd;

  boost::progress_display* progressBar = 0;
  if ( progress ) {
    Repository::clog() << "exploring " << process();
    progressBar = new boost::progress_display(theExplorationSteps,Repository::clog());
  }
  std::set<SimpleCellGrid*> newCells;
  
  if ( pre_adaption_splits().empty() &&
       (theLuminositySplits || theChannelSplits || theAllChannelSplits) ) {
    const StandardEventHandler& eh = *eventHandler();
    const StandardXComb& xc = *eh.xCombs()[bin()];
    the_pre_adaption_splits.resize(dimension(),0);
    const pair<int,int>& pdims = xc.partonDimensions();
    if ( theLuminositySplits && dimension() >= pdims.first + pdims.second ) {
      for ( int n = 0; n < pdims.first; ++n )
	the_pre_adaption_splits[n] = theLuminositySplits;
      for ( int n = dimension() - pdims.second; n < dimension(); ++n )
	the_pre_adaption_splits[n] = theLuminositySplits;
    }
    if ( theChannelSplits && xc.diagrams().size() &&
	 dimension() > pdims.first + pdims.second ) {
      the_pre_adaption_splits[pdims.first] = theChannelSplits;
    }
    if ( theAllChannelSplits && xc.diagrams().size() > 1 &&
	 dimension() > pdims.first + pdims.second ) {
      the_pre_adaption_splits[pdims.first] = xc.diagrams().size() - 1;
    }
  }
  
  for(int splitdim=0; splitdim<min(dimension(),(int)pre_adaption_splits().size());splitdim++)
      SimpleCellGrid::splitter(splitdim,pre_adaption_splits()[splitdim]);
  
  SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells,Repository::clog());

  
  bool notAll = false;
  for ( std::size_t step = 1; step < theExplorationSteps; ++step ) {
    newCells.clear();
    SimpleCellGrid::adapt(theGain,theEpsilon,newCells);
    if ( progressBar )
      ++(*progressBar);
    if ( newCells.empty() ) {
      notAll = true;
      break;
    }
    SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells,Repository::clog());
  }

  if ( progressBar )
    ++(*progressBar);

  SimpleCellGrid::setWeights();
  SimpleCellGrid::updateIntegral();
  SimpleCellGrid::minimumSelection(theMinimumSelection);

  if ( progressBar ) {
    if ( notAll )
      cout << "\n" << flush;
    delete progressBar;
  }

  unsigned long points = initialPoints();
  for ( unsigned long k = 0; k < nIterations(); ++k ) {
    runIteration(points,progress);
    if ( k < nIterations() - 1 ) {
      points = (unsigned long)(points*enhancementFactor());
      adapt();
      nextIteration();
    }
  }

  didReadGrids();
  isInitialized();

}

void CellGridSampler::finalize(bool) {
  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);
  if (randomNumberString()!="")  
  for ( map<RandomNumberIndex,pair<RandomNumberHistogram,double> >::
    const_iterator b = RandomNumberHistograms.begin();
    b != RandomNumberHistograms.end(); ++b ) {
    b->second.first.dump(randomNumberString(), b->first.first,shortprocess(),b->first.second);
  }


}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void CellGridSampler::persistentOutput(PersistentOStream & os) const {
  os << theExplorationPoints << theExplorationSteps
     << theGain << theEpsilon << theMinimumSelection
     << the_pre_adaption_splits
     << theLuminositySplits << theChannelSplits
     << theAllChannelSplits << theUnweightCells;
}

void CellGridSampler::persistentInput(PersistentIStream & is, int) {
  is >> theExplorationPoints >> theExplorationSteps
     >> theGain >> theEpsilon >> theMinimumSelection
     >> the_pre_adaption_splits
     >> theLuminositySplits >> theChannelSplits
     >> theAllChannelSplits >> theUnweightCells;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<CellGridSampler,BinSampler>
  describeHerwigCellGridSampler("Herwig::CellGridSampler", "HwSampling.so");

void CellGridSampler::Init() {

  static ClassDocumentation<CellGridSampler> documentation
    ("CellGridSampler samples XCombs bins using CellGrids.");

  static Parameter<CellGridSampler,size_t> interfaceExplorationPoints
    ("ExplorationPoints",
     "The number of points to use for cell exploration.",
     &CellGridSampler::theExplorationPoints, 1000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler,size_t> interfaceExplorationSteps
    ("ExplorationSteps",
     "The number of exploration steps to perform.",
     &CellGridSampler::theExplorationSteps, 8, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler,double> interfaceGain
    ("Gain",
     "The gain factor used for adaption.",
     &CellGridSampler::theGain, 0.3, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<CellGridSampler,double> interfaceEpsilon
    ("Epsilon",
     "The efficieny threshold used for adaption.",
     &CellGridSampler::theEpsilon, 0.01, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<CellGridSampler,double> interfaceMinimumSelection
    ("MinimumSelection",
     "The minimum cell selection probability.",
     &CellGridSampler::theMinimumSelection, 0.0001, 0.0, 1.0,
     false, false, Interface::limited);
    
  static ParVector<CellGridSampler,int> interfacethe_pre_adaption_splits
    ("preadaptionsplit",
     "The splittings for each dimension befor adaption.",
     &CellGridSampler::the_pre_adaption_splits, 1., -1, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler,int> interfaceLuminositySplits
    ("LuminositySplits",
     "",
     &CellGridSampler::theLuminositySplits, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler,int> interfaceChannelSplits
    ("ChannelSplits",
     "",
     &CellGridSampler::theChannelSplits, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<CellGridSampler,bool> interfaceAllChannelSplits
    ("AllChannelSplits",
     "",
     &CellGridSampler::theAllChannelSplits, false, false, false);
  static SwitchOption interfaceAllChannelSplitsYes
    (interfaceAllChannelSplits,
     "Yes",
     "",
     true);
  static SwitchOption interfaceAllChannelSplitsNo
    (interfaceAllChannelSplits,
     "No",
     "",
     false);

  static Switch<CellGridSampler,bool> interfaceUnweightCells
    ("UnweightCells",
     "",
     &CellGridSampler::theUnweightCells, true, false, false);
  static SwitchOption interfaceUnweightCellsYes
    (interfaceUnweightCells,
     "Yes",
     "",
     true);
  static SwitchOption interfaceUnweightCellsNo
    (interfaceUnweightCells,
     "No",
     "",
     false);


}

