// -*- C++ -*-
//
// SimpleCellGridSampler.cpp is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleCellGridSampler class.
//

#include "SimpleCellGridSampler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <boost/progress.hpp>

#include "CellGridSampler.h"

using namespace Herwig;
using namespace ExSample;

SimpleCellGridSampler::SimpleCellGridSampler() 
  : BinSampler(), SimpleCellGrid(),
    theExplorationPoints(5000), theExplorationSteps(10),
    theNIterations(1), theEnhancementFactor(2.0),
    theGain(0.5), theMinimumSelection(0.1),
    theLastNPoints(0) {}

SimpleCellGridSampler::~SimpleCellGridSampler() {}

IBPtr SimpleCellGridSampler::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleCellGridSampler::fullclone() const {
  return new_ptr(*this);
}

void SimpleCellGridSampler::generate(bool noMaxInfo) {

  UseRandom rnd;
  pair<double,double> weights = SimpleCellGrid::generate(rnd,*this,lastPoint());

  assert(weights.second != 0.0);
  double w = SimpleCellGrid::integral()*weights.first/weights.second;

  if ( abs(w) > maxWeight() && !noMaxInfo ) {
    double old = maxWeight();
    select(w);
    throw NewMaximum(old,abs(w));
  }

  select(w);

}

void SimpleCellGridSampler::adapt() {
  UseRandom rnd;
  set<SimpleCellGrid*> newCells;
  SimpleCellGrid::adapt(theGain,newCells);
  SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells);
  SimpleCellGrid::updateIntegral();
  SimpleCellGrid::minimumSelection(theMinimumSelection);
}

void SimpleCellGridSampler::initialize(bool progress) {

  bool haveGrid = false;
  Ptr<CellGridSampler>::tptr sampler = 
    dynamic_ptr_cast<Ptr<CellGridSampler>::tptr>(eventHandler()->sampler());
  if ( !sampler )
    throw Exception() << "SimpleCellGridSampler: Need to be contained in a CellGridSampler object."
		      << Exception::abortnow;

  list<XML::Element>::iterator git = sampler->grids().children().begin();
  for ( ; git != sampler->grids().children().end(); ++git ) {
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
    sampler->grids().erase(git);
  }

  lastPoint().resize(dimension());

  if ( initialized() ) {
    if ( !haveGrid )
      throw Exception() << "SimpleCellGridSampler: Require existing grid when starting to run."
			<< Exception::abortnow;
    return;
  }

  if ( haveGrid ) {
    runIteration(initialPoints(),progress);
    isInitialized();
    XML::Element grid = SimpleCellGrid::toXML();
    grid.appendAttribute("process",id());
    sampler->grids().append(grid);
    return;
  }

  SimpleCellGrid::boundaries(vector<double>(dimension(),0.0),vector<double>(dimension(),1.0));
  SimpleCellGrid::weightInformation().resize(dimension());

  UseRandom rnd;
  boost::progress_display* progressBar = 0;
  if ( progress ) {
    cout << "exploring " << process() << "\n"
	 << "(id " << id() << ")";
    progressBar = new boost::progress_display(theExplorationSteps,cout);
  }

  std::set<SimpleCellGrid*> newCells;
  bool notAll = false;
  for ( std::size_t step = 0; step < theExplorationSteps; ++step ) {
    SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells);
    if ( progressBar )
      ++(*progressBar);
    newCells.clear();
    SimpleCellGrid::adapt(theGain,newCells);
    if ( newCells.empty() ) {
      notAll = true;
      break;
    }
  }
  SimpleCellGrid::updateIntegral();
  SimpleCellGrid::minimumSelection(theMinimumSelection);

  if ( progressBar ) {
    if ( notAll )
      cout << "\n" << flush;
    delete progressBar;
  }

  theLastNPoints = initialPoints();
  runIteration(theLastNPoints,progress);
  theLastNPoints = (unsigned long)(theLastNPoints*theEnhancementFactor);
  if ( theNIterations == 1 ) {
    isInitialized();
    XML::Element grid = SimpleCellGrid::toXML();
    grid.appendAttribute("process",id());
    sampler->grids().append(grid);
    return;
  }
  adapt();
  nextIteration();
  for ( unsigned long k = 1; k < theNIterations-1; ++k ) {
    runIteration(theLastNPoints,progress);
    theLastNPoints = (unsigned long)(theLastNPoints*theEnhancementFactor);
    adapt();
    nextIteration();
  }
  runIteration(theLastNPoints,progress);
  theLastNPoints = (unsigned long)(theLastNPoints*theEnhancementFactor);
  adapt();
  isInitialized();

  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler->grids().append(grid);

}

void SimpleCellGridSampler::finalize(bool) {

  Ptr<CellGridSampler>::tptr sampler = 
    dynamic_ptr_cast<Ptr<CellGridSampler>::tptr>(eventHandler()->sampler());
  if ( !sampler )
    throw Exception() << "SimpleCellGridSampler: Need to be contained in a CellGridSampler object."
		      << Exception::abortnow;

  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler->grids().append(grid);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void SimpleCellGridSampler::persistentOutput(PersistentOStream & os) const {
  os << theExplorationPoints << theExplorationSteps
     << theNIterations << theEnhancementFactor
     << theGain << theMinimumSelection
     << theLastNPoints;
}

void SimpleCellGridSampler::persistentInput(PersistentIStream & is, int) {
  is >> theExplorationPoints >> theExplorationSteps
     >> theNIterations >> theEnhancementFactor
     >> theGain >> theMinimumSelection
     >> theLastNPoints;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SimpleCellGridSampler,BinSampler>
  describeHerwigSimpleCellGridSampler("Herwig::SimpleCellGridSampler", "HwExsample2.so");

void SimpleCellGridSampler::Init() {

  static ClassDocumentation<SimpleCellGridSampler> documentation
    ("SimpleCellGridSampler samples XCombs bins using CellGrids.");

  static Parameter<SimpleCellGridSampler,size_t> interfaceExplorationPoints
    ("ExplorationPoints",
     "The number of points to use for cell exploration.",
     &SimpleCellGridSampler::theExplorationPoints, 5000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<SimpleCellGridSampler,size_t> interfaceExplorationSteps
    ("ExplorationSteps",
     "The number of exploration steps to perform.",
     &SimpleCellGridSampler::theExplorationSteps, 10, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<SimpleCellGridSampler,size_t> interfaceNIterations
    ("NIterations",
     "The number of iterations to perform initially.",
     &SimpleCellGridSampler::theNIterations, 1, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<SimpleCellGridSampler,double> interfaceEnhancementFactor
    ("EnhancementFactor",
     "The enhancement factor for the number of points in the next iteration.",
     &SimpleCellGridSampler::theEnhancementFactor, 2.0, 1.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<SimpleCellGridSampler,double> interfaceGain
    ("Gain",
     "The gain factor used for adaption.",
     &SimpleCellGridSampler::theGain, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<SimpleCellGridSampler,double> interfaceMinimumSelection
    ("MinimumSelection",
     "The minimum cell selection probability.",
     &SimpleCellGridSampler::theMinimumSelection, 0.1, 0.0, 1.0,
     false, false, Interface::limited);


}

