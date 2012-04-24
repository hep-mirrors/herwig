// -*- C++ -*-
//
// ProjectingSampler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProjectingSampler class.
//

#include "ProjectingSampler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ProjectingSampler::ProjectingSampler() 
  : theFirstIteration(true), theNIterations(4),
    theEnhancementFactor(2.0), theNBins(8),
    theEpsilon(0.5), theLastNPoints(0),
    theWeightThreshold(0.001) {}

ProjectingSampler::~ProjectingSampler() {}

IBPtr ProjectingSampler::clone() const {
  return new_ptr(*this);
}

IBPtr ProjectingSampler::fullclone() const {
  return new_ptr(*this);
}

void ProjectingSampler::select(double weight) {
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    if ( theFirstIteration ) {
      theProjections[k].bin(lastPoint()[k]);
    }
    theProjections[k].select(theLastValue);
  }
  GeneralStatistics::select(weight);
}

void ProjectingSampler::accept() {
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    theProjections[k].accept();
  }
  GeneralStatistics::accept();
}

void ProjectingSampler::reject() {
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    theProjections[k].reject();
  }
  GeneralStatistics::reject();
}

void ProjectingSampler::generate(bool noMaxInfo) {
  double w = 1.;
  if ( !theFirstIteration ) {
    for ( size_t k = 0; k < lastPoint().size(); ++k ) {
      w *= theProjections[k].sample(lastPoint()[k]);
    }
  } else {
    for ( size_t k = 0; k < lastPoint().size(); ++k )
      lastPoint()[k] = UseRandom::rnd();
  }
  try {
    theLastValue = theEventHandler->dSigDR(lastPoint()) / nanobarn;
    w *= theLastValue;
  } catch (Veto&) {
    theLastValue = 0.0;
    w = 0.0;
  } catch (...) {
    throw;
  }
  if ( abs(w) > maxWeight() && !noMaxInfo ) {
    double old = maxWeight();
    select(w);
    throw NewMaximum(old,abs(w));
  }
  select(w);
  if ( !noMaxInfo &&
       selectedPoints() == theLastNPoints ) {
    theLastNPoints = (unsigned long)(theLastNPoints*theEnhancementFactor);
    adapt();
    throw UpdateCrossSections();
  }
}

void ProjectingSampler::initialize(bool progress) {
  if ( initialized() )
    return;
  lastPoint().resize(dimension());
  theProjections.resize(dimension(),BinnedStatistics(theNBins,theWeightThreshold));
  theLastNPoints = initialPoints();
  theFirstIteration = true;
  runIteration(theLastNPoints,progress);
  theLastNPoints = (unsigned long)(theLastNPoints*theEnhancementFactor);
  adapt();
  if ( theNIterations == 1 ) {
    isInitialized();
    return;
  }
  nextIteration();
  theFirstIteration = false;
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
}

struct ProjectingAdaptor {

  double variance;
  double epsilon;

  double importanceMeasure(const GeneralStatistics& s) const {
    return s.averageAbsWeight();
  }

  bool adapt(const GeneralStatistics& s) const {
    return s.averageAbsWeightVariance()/variance > epsilon;
  }

};

void ProjectingSampler::adapt() {

  ProjectingAdaptor adaptor;
  adaptor.variance = 0.;
  for ( vector<BinnedStatistics>::iterator s =
	  theProjections.begin(); s != theProjections.end(); ++s ) {
    double variance = 0.;
    for ( map<double,GeneralStatistics>::const_iterator k = 
	    s->statistics().begin(); k != s->statistics().end(); ++k ) {
      variance += k->second.averageAbsWeightVariance();
    }
    adaptor.variance += variance/s->statistics().size();
  }
  adaptor.variance /= theProjections.size();
  adaptor.epsilon = theEpsilon;

  size_t count = 0;
  for ( vector<BinnedStatistics>::iterator s =
	  theProjections.begin(); s != theProjections.end(); ++s, ++count ) {
    s->adapt(adaptor);
  }

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ProjectingSampler::persistentOutput(PersistentOStream & os) const {
  os << theFirstIteration << theNIterations << theEnhancementFactor 
     << theNBins << theEpsilon << theLastNPoints
     << theProjections << theWeightThreshold;
}

void ProjectingSampler::persistentInput(PersistentIStream & is, int) {
  is >> theFirstIteration >> theNIterations >> theEnhancementFactor 
     >> theNBins >> theEpsilon >> theLastNPoints
     >> theProjections >> theWeightThreshold;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ProjectingSampler,Herwig::BinSampler>
  describeHerwigProjectingSampler("Herwig::ProjectingSampler", "HwExsample2.so");

void ProjectingSampler::Init() {

  static ClassDocumentation<ProjectingSampler> documentation
    ("ProjectingSampler does adaption from projections of the integrand.");


  static Parameter<ProjectingSampler,unsigned long> interfaceNIterations
    ("NIterations",
     "The number of iterations to perform initially.",
     &ProjectingSampler::theNIterations, 4, 1, 0,
     false, false, Interface::lowerlim);


  static Parameter<ProjectingSampler,double> interfaceEnhancementFactor
    ("EnhancementFactor",
     "The enhancement factor for the number of points in the next iteration.",
     &ProjectingSampler::theEnhancementFactor, 2.0, 1.0, 0,
     false, false, Interface::lowerlim);


  static Parameter<ProjectingSampler,unsigned int> interfaceNBins
    ("NBins",
     "The number of projection bins to consider initially.",
     &ProjectingSampler::theNBins, 8, 1, 0,
     false, false, Interface::lowerlim);


  static Parameter<ProjectingSampler,double> interfaceEpsilon
    ("Epsilon",
     "The adaption threshold.",
     &ProjectingSampler::theEpsilon, 0.5, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ProjectingSampler,double> interfaceWeightThreshold
    ("WeightThreshold",
     "The minimum weight per bin in units of the average weight.",
     &ProjectingSampler::theWeightThreshold, 0.001, 0.0, 0,
     false, false, Interface::lowerlim);

}

