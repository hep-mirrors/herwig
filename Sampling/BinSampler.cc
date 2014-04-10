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
#include "ThePEG/Repository/Repository.h"

#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <boost/progress.hpp>

#include "GeneralSampler.h"

using namespace Herwig;

BinSampler::BinSampler() 
  : MultiIterationStatistics(), 
    theBias(1.),
    theWeighted(false),
    theInitialPoints(1000000),
    theNIterations(1),
    theEnhancementFactor(1.0),
    theReferenceWeight(1.0),
    theBin(-1),
    theInitialized(false) {}

BinSampler::~BinSampler() {}

IBPtr BinSampler::clone() const {
  return new_ptr(*this);
}

IBPtr BinSampler::fullclone() const {
  return new_ptr(*this);
}

void BinSampler::sampler(Ptr<GeneralSampler>::tptr s) {
  theSampler = s;
}

Ptr<GeneralSampler>::tptr BinSampler::sampler() const {
  return theSampler;
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

string BinSampler::shortprocess() const {
  ostringstream os("");
  const StandardEventHandler& eh = *theEventHandler;
  const StandardXComb& xc = *eh.xCombs()[theBin];
  os << xc.mePartonData()[0]->id() << " "
     << xc.mePartonData()[1]->id() << " : ";
  for ( cPDVector::const_iterator pid =
	  xc.mePartonData().begin() + 2;
	pid != xc.mePartonData().end(); ++pid )
    os << (**pid).id() << " ";
  return os.str();
}

string BinSampler::id() const {
  ostringstream os("");
  const StandardEventHandler& eh = *theEventHandler;
  const StandardXComb& xc = *eh.xCombs()[theBin];
  string name = xc.matrixElement()->name();
  string::size_type i = name.find_first_of("[");
  string nameFirst = name.substr(0,i);
  i = name.find_first_of("]");
  string nameSecond = name.substr(i+1);
  os << nameFirst << nameSecond << ":";
  for ( cPDVector::const_iterator pid =
	  xc.mePartonData().begin();
	pid != xc.mePartonData().end(); ++pid )
    os << (**pid).id() << (pid != (--xc.mePartonData().end()) ? "," : "");
  return os.str();
}

double BinSampler::generate() {
  double w = 1.;
//  cout<<"\npoint: ";
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    lastPoint()[k] = UseRandom::rnd();
//    cout<<lastPoint()[k]<<" ";
  }
  try {
    w = eventHandler()->dSigDR(lastPoint()) / nanobarn;
  } catch (Veto&) {
    w = 0.0;
  } catch (...) {
    throw;
  }
  if (randomNumberString()!="") 
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    RandomNumberHistograms[RandomNumberIndex(id(),k)].first.book(lastPoint()[k],w);
    RandomNumberHistograms[RandomNumberIndex(id(),k)].second+=w;
  }

  
  
  if ( !weighted() && initialized() ) {
    double p = min(abs(w),referenceWeight())/referenceWeight();
    double sign = w >= 0. ? 1. : -1.;
    if ( p < 1 && UseRandom::rnd() > p )
      w = 0.;
    else
      w = sign*max(abs(w),referenceWeight());
  }
  select(w);
  if ( w != 0.0 )
    accept();
  return w;
}

void BinSampler::runIteration(unsigned long points, bool progress) {

  boost::progress_display* progressBar = 0;
  if ( progress ) {
    Repository::clog() << "integrating " << process() << " , iteration "
		       << (iterations().size() + 1);
    progressBar = new boost::progress_display(points,Repository::clog());
  }

  for ( unsigned long k = 0; k < points; ++k ) {

    generate();

    if ( progress ) {
      ++(*progressBar);
    }

  }

  if ( progress ) {
    Repository::clog() << "integrated ( " 
		       << averageWeight() << " +/- " << sqrt(averageWeightVariance())
		       << " ) nb\nepsilon = "
		       << (abs(maxWeight()) != 0. ? averageAbsWeight()/abs(maxWeight()) : 0.);
    if ( !iterations().empty() )
      Repository::clog() << " chi2 = " << chi2();
    Repository::clog() << "\n";
    Repository::clog() << "--------------------------------------------------------------------------------\n";
  }

  if ( progressBar )
    delete progressBar;

}

void BinSampler::initialize(bool progress) {
  lastPoint().resize(dimension());
  if (randomNumberString()!="") 
  for(size_t i=0;i<lastPoint().size();i++){
     RandomNumberHistograms[RandomNumberIndex(id(),i)] = make_pair( RandomNumberHistogram(),0.);
  }
  if ( initialized() )
    return;
  if ( !sampler()->grids().children().empty() ) {
    nIterations(1);
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
  isInitialized();
}


void BinSampler::finalize(bool){
  if (theRandomNumbers!="")  
    for ( map<RandomNumberIndex,pair<RandomNumberHistogram,double> >::
	  const_iterator b = RandomNumberHistograms.begin();
	b != RandomNumberHistograms.end(); ++b ) {
       b->second.first.dump(randomNumberString(), b->first.first,shortprocess(),b->first.second);
  }

}




BinSampler::RandomNumberHistogram::
RandomNumberHistogram(double low, 
		     double up, 
		     unsigned int nbins) 
  : lower(low) {
  nbins = nbins + 1;

  double c = up / (nbins-1.);

  for ( unsigned int k = 1; k < nbins; ++k ) {
    bins[low+c*k] = 0.;
    binsw1[low+c*k] = 0.;
  }

}


void BinSampler::RandomNumberHistogram::
dump(const std::string& folder,const std::string& prefix, const std::string& process, 
     const int NR) const {
  ostringstream fname("");
  std::string prefix2;
  std::string prefix3=prefix;
  std::remove_copy(prefix.begin(), prefix.end(), std::back_inserter(prefix2), '.');
  prefix3=prefix2;prefix2.clear();
  std::remove_copy(prefix3.begin(), prefix3.end(), std::back_inserter(prefix2), ':');
    prefix3=prefix2;prefix2.clear();
  std::remove_copy(prefix3.begin(), prefix3.end(), std::back_inserter(prefix2), ',');
  fname << "RN-"<< NR ;
  ofstream out((folder+"/"+prefix2+fname.str()+".dat").c_str());
  double sumofweights=0.;
  for ( map<double,double >::const_iterator b = bins.begin();b != bins.end(); ++b )
       sumofweights+=b->second;  
  for ( map<double,double >::const_iterator b = bins.begin();
	b != bins.end(); ++b ) {
      out << " " << b->first
	  << " " << b->second/sumofweights*100.
	  << "\n" << flush;
  } 
  sumofweights=0.;
  for ( map<double,double >::const_iterator b = binsw1.begin();b != binsw1.end(); ++b )
       sumofweights+=b->second;  
  
  ofstream out2((folder+"/"+prefix2+fname.str()+"-w=1.dat").c_str());
  for ( map<double,double >::const_iterator b = binsw1.begin();
	b != binsw1.end(); ++b ) {
      out2 << " " << b->first
	  << " " << b->second/sumofweights*100.
	  << "\n" << flush;
  }
  double xmin = -0.01;
  double xmax = 1.01;
  ofstream gpout((folder+"/"+prefix2+fname.str()+".gp").c_str());
  gpout << "set terminal epslatex color solid\n"
      << "set output '" << prefix2+fname.str() << "-plot.tex'\n"
      << "set xrange [" << xmin << ":" << xmax << "]\n";
    gpout << "set xlabel 'rn "<<NR <<"' \n";
    gpout << "set size 0.5,0.6\n";
    gpout << "plot '" << prefix2+fname.str()
    << ".dat' u ($1):($2)  w boxes  lc rgbcolor \"blue\" t '{\\tiny "<<process <<"}',";
    gpout << " '" << prefix2+fname.str();
    gpout << "-w=1.dat' u ($1):($2)  w boxes  lc rgbcolor \"red\" t '{\\tiny "<<process <<":w=1}';";
  gpout << "reset\n";
}








// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void BinSampler::persistentOutput(PersistentOStream & os) const {
  MultiIterationStatistics::put(os);
  os << theBias << theWeighted << theInitialPoints << theNIterations 
     << theEnhancementFactor << theReferenceWeight
     << theBin << theInitialized << theLastPoint
     << theEventHandler << theSampler << theRandomNumbers;
}

void BinSampler::persistentInput(PersistentIStream & is, int) {
  MultiIterationStatistics::get(is);
  is >> theBias >> theWeighted >> theInitialPoints >> theNIterations 
     >> theEnhancementFactor >> theReferenceWeight
     >> theBin >> theInitialized >> theLastPoint
     >> theEventHandler >> theSampler >> theRandomNumbers;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<BinSampler,MultiIterationStatistics>
  describeHerwigBinSampler("Herwig::BinSampler", "HwSampling.so");

void BinSampler::Init() {

  static ClassDocumentation<BinSampler> documentation
    ("BinSampler samples XCombs bins. This default implementation performs flat MC integration.");

  static Parameter<BinSampler,unsigned long> interfaceInitialPoints
    ("InitialPoints",
     "The number of points to use for initial integration.",
     &BinSampler::theInitialPoints, 1000000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<BinSampler,size_t> interfaceNIterations
    ("NIterations",
     "The number of iterations to perform initially.",
     &BinSampler::theNIterations, 1, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<BinSampler,double> interfaceEnhancementFactor
    ("EnhancementFactor",
     "The enhancement factor for the number of points in the next iteration.",
     &BinSampler::theEnhancementFactor, 2.0, 1.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<BinSampler,string> interfaceRandomNumbers
    ("RandomNumbers",
     "Prefix for distributions of the random numbers.",
     &BinSampler::theRandomNumbers, "",
     false, false);
}

