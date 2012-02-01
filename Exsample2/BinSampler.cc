// -*- C++ -*-
//
// BinSampler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BinSampler class.
//

#include "BinSampler.h"
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

using namespace Herwig;

BinSampler::BinSampler() 
  : Interfaced(), MultiIterationStatistics(), 
    theInitialPoints(1000000), theBin(-1),
    theInitialized(false) {}

BinSampler::~BinSampler() {}

IBPtr BinSampler::clone() const {
  return new_ptr(*this);
}

IBPtr BinSampler::fullclone() const {
  return new_ptr(*this);
}

string BinSampler::process() const {
  ostringstream os("");
  const StandardEventHandler& eh = *theEventHandler;
  const StandardXComb& xc = *eh.xCombs()[theBin];
  os << xc.matrixElement()->name() << " : ";
  os << xc.mePartonData()[0]->PDGName() << " "
     << xc.mePartonData()[1]->PDGName() << " -> ";
  for ( cPDVector::const_iterator pid =
	  xc.mePartonData().begin() + 2;
	pid != xc.mePartonData().end(); ++pid )
    os << (**pid).PDGName() << " ";
  return os.str();
}

void BinSampler::generate(bool noMaxInfo) {
  double w = 1.;
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    lastPoint()[k] = UseRandom::rnd();
  }
  try {
    w = theEventHandler->dSigDR(lastPoint()) / nanobarn;
  } catch (Veto&) {
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
}

void BinSampler::runIteration(unsigned long points, bool progress) {

  boost::progress_display* progressBar = 0;
  if ( progress ) {
    cout << "integrating " << process() << " , iteration "
	 << (iterations().size() + 1);
    progressBar = new boost::progress_display(points,cout);
  }

  for ( unsigned long k = 0; k < points; ++k ) {

    generate(true);

    if ( progress ) {
      ++(*progressBar);
    }

  }

  if ( !iterations().empty() )
    chi2(iterations().back());

  if ( progress ) {
    cout << "integrated ( " 
	 << averageWeight() << " +/- " << sqrt(averageWeightVariance())
	 << " ) nb\nepsilon = "
	 << averageAbsWeight()/abs(maxWeight());
    if ( chi2() >= 0. )
      cout << " chi2 = " << chi2();
    cout << "\n";
    cout << "--------------------------------------------------------------------------------\n";
  }

  if ( progressBar )
    delete progressBar;

}

void BinSampler::initialize(bool progress) {
  if ( initialized() )
    return;
  lastPoint().resize(dimension());
  runIteration(initialPoints(),progress);
  isInitialized();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void BinSampler::persistentOutput(PersistentOStream & os) const {
  MultiIterationStatistics::put(os);
  os << theInitialPoints << theBin << theInitialized << theLastPoint;
}

void BinSampler::persistentInput(PersistentIStream & is, int) {
  MultiIterationStatistics::get(is);
  is >> theInitialPoints >> theBin >> theInitialized >> theLastPoint;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<BinSampler,Interfaced>
  describeHerwigBinSampler("Herwig::BinSampler", "HwExsample2.so");

void BinSampler::Init() {

  static ClassDocumentation<BinSampler> documentation
    ("BinSampler samples XCombs bins. This default implementation performs flat MC integration.");


  static Parameter<BinSampler,unsigned long> interfaceInitialPoints
    ("InitialPoints",
     "The number of points to use for initial integration.",
     &BinSampler::theInitialPoints, 1000000, 1, 0,
     false, false, Interface::lowerlim);

}

