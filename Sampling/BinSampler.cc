// -*- C++ -*-
//
// BinSampler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

#include "ThePEG/Utilities/ColourOutput.h"

#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include "Herwig/Utilities/Progress.h"

#include "GeneralSampler.h"

using namespace Herwig;

BinSampler::BinSampler() 
  : MultiIterationStatistics(), 
    theBias(1.),
    theWeighted(false),
    theInitialPoints(1000000),
    theNIterations(1),
    theEnhancementFactor(1.0),
    theNonZeroInPresampling(false),
    theHalfPoints(false),
    theMaxNewMax(30),
    theReferenceWeight(1.0),
    theBin(-1),
    theInitialized(false),
    theRemapperPoints(0),
    theRemapChannelDimension(false),
    theLuminosityMapperBins(0),
    theGeneralMapperBins(0),
    theRemapperMinSelection(0.00001),
    theIntegrated(false),
    theRemappersFilled(false),
    theHasGrids(false),
    theKappa(1.){}

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

double BinSampler::evaluate(vector<double> p,
			    bool remap) {
  double w = 1.0;
  if ( remap && !remappers.empty() ) {
    for ( size_t k = 0; k < p.size(); ++k ) {
      map<size_t,Remapper>::const_iterator r =
	remappers.find(k);
      if ( r != remappers.end() ) {
	pair<double,double> f = r->second.generate(p[k]);
	p[k] = f.first;
	w /= f.second;
      }
    }
  }
  try {
    w *= eventHandler()->dSigDR(p) / nanobarn;
  } catch (Veto&) {
    w = 0.0;
  } catch (...) {
    throw;
  }
  if (randomNumberString()!="") 
    for ( size_t k = 0; k < p.size(); ++k ) {
      RandomNumberHistograms[RandomNumberIndex(id(),k)].first.book(p[k],w);
      RandomNumberHistograms[RandomNumberIndex(id(),k)].second+=w;
    }
  return w;
}

double BinSampler::generate() {
  double w = 1.;
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    lastPoint()[k] = UseRandom::rnd();
  }
  try {
    w = evaluate(lastPoint());
  } catch (Veto&) {
    w = 0.0;
  } catch (...) {
    throw;
  }
 
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

void BinSampler::fillRemappers(bool progress) {

  if ( remappers.empty() )
    return;

  unsigned long nanPoints = 0;

  progress_display* progressBar = nullptr;
  if ( progress ) {
    Repository::clog() << "warming up " << ANSI::red << process() << ANSI::reset;
    progressBar = new progress_display{ theRemapperPoints, Repository::clog() };
  }

  unsigned long countzero =0;
  for ( unsigned long k = 0; k < theRemapperPoints; ++k,++countzero ) {
    if (countzero>=theRemapperPoints)break;
    double w = 1.;
    for ( size_t j = 0; j < lastPoint().size(); ++j ) {
      lastPoint()[j] = UseRandom::rnd();
    }
    try {
      w = evaluate(lastPoint(),false);
    } catch (Veto&) {
      w = 0.0;
    } catch (...) {
      throw;
    }

    if ( ! isfinite(w) )
      ++nanPoints;
    
    if ( theNonZeroInPresampling &&  w==0. ){
      k--;
      continue; 
    }

    if ( w != 0.0 ) {
      countzero=0;
      for ( map<size_t,Remapper>::iterator r = remappers.begin();
	    r != remappers.end(); ++r )
	r->second.fill(lastPoint()[r->first],w);
    }

    if ( progressBar )
      ++(*progressBar);

  }

  if ( progressBar ) {
    delete progressBar;
  }

  if ( nanPoints ) {
    Repository::clog() << "Warning: " << nanPoints 
		       << " out of " << theRemapperPoints << " points with nan or inf "
		       << "weight encountered while filling remappers.\n" << flush;
  }

}

void BinSampler::saveIntegrationData() const {

  XML::Element stats = MultiIterationStatistics::toXML();
  stats.appendAttribute("process",id());

  sampler()->grids().append(stats);

}

void BinSampler::readIntegrationData() {

  if ( theIntegrated )
    return;

  bool haveStats = false;

  list<XML::Element>::iterator sit = sampler()->grids().children().begin();
  for ( ; sit != sampler()->grids().children().end(); ++sit ) {
    if ( sit->type() != XML::ElementTypes::Element )
      continue;
    if ( sit->name() != "MultiIterationStatistics" )
      continue;
    string proc;
    sit->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveStats = true;
      break;
    }
  }

  if ( haveStats ) {
    MultiIterationStatistics::fromXML(*sit);
    sampler()->grids().erase(sit);
    theIntegrated = true;
  } else {
    throw Exception()
      << "\n---------------------------------------------------\n\n"
      << "Expected integration data.\n\n"
      << "* When using the build setup make sure the integrate command has been run.\n\n"
      << "* Check the [EventGenerator].log file for further information.\n\n"
      << "* Make sure that the Herwig folder can be found and that it contains a HerwigGrids.xml file.\n\n"
      << "* If you have split the integration jobs, make sure that each integration job was finished.\n"
      << "  Afterwards delete the global HerwigGrids.xml file in the Herwig subfolder\n"
      << "  to automatically create an updated version of the global HerwigGrids.xml file.\n\n"
      << "---------------------------------------------------\n"
      << Exception::abortnow;
  }

}

void BinSampler::saveRemappers() const {

  if ( remappers.empty() )
    return;

  XML::Element maps(XML::ElementTypes::Element,"Remappers");
  maps.appendAttribute("process",id());

  for ( map<size_t,Remapper>::const_iterator r = remappers.begin();
	r != remappers.end(); ++r ) {
    XML::Element rmap = r->second.toXML();
    rmap.appendAttribute("dimension",r->first);
    maps.append(rmap);
  }

  sampler()->grids().append(maps);

}

void BinSampler::setupRemappers(bool progress) {

  if ( !theRemapperPoints )
    return;

  if ( theRemappersFilled )
    return;

  lastPoint().resize(dimension());

  bool haveGrid = false;

  list<XML::Element>::iterator git = sampler()->grids().children().begin();
  for ( ; git != sampler()->grids().children().end(); ++git ) {
    if ( git->type() != XML::ElementTypes::Element )
      continue;
    if ( git->name() != "Remappers" )
      continue;
    string proc;
    git->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveGrid = true;
      break;
    }
  }

  if ( haveGrid ) {
    for ( list<XML::Element>::iterator cit = git->children().begin();
	  cit != git->children().end(); ++cit ) {
      if ( cit->type() != XML::ElementTypes::Element )
	continue;
      if ( cit->name() != "Remapper" )
	continue;
      size_t dimension = 0;
      cit->getFromAttribute("dimension",dimension);
      remappers[dimension].fromXML(*cit);
    }
    sampler()->grids().erase(git);
  }

  if ( !haveGrid ) {

    const StandardEventHandler& eh = *eventHandler();
    const StandardXComb& xc = *eh.xCombs()[bin()];
    const pair<int,int>& pdims = xc.partonDimensions();

    set<int> remapped;

    if ( theRemapChannelDimension && xc.diagrams().size() > 1 &&
	 dimension() > pdims.first + pdims.second ) {
      remappers[pdims.first] = Remapper(xc.diagrams().size(),theRemapperMinSelection,false);
      remapped.insert(pdims.first);
    }

    if ( theLuminosityMapperBins > 1 && dimension() >= pdims.first + pdims.second ) {
      for ( int n = 0; n < pdims.first; ++n ) {
	remappers[n] = Remapper(theLuminosityMapperBins,theRemapperMinSelection,true);
	remapped.insert(n);
      }
      for ( int n = dimension() - pdims.second; n < dimension(); ++n ) {
	remappers[n] = Remapper(theLuminosityMapperBins,theRemapperMinSelection,true);
	remapped.insert(n);
      }
    }

    if ( theGeneralMapperBins > 1 ) {
      for ( int n = 0; n < dimension(); n++ ) {
	if ( remapped.find(n) == remapped.end() ) {
	  remappers[n] = Remapper(theGeneralMapperBins,theRemapperMinSelection,true);
	  remapped.insert(n);
	}
      }
    }

    fillRemappers(progress);

    for ( map<size_t,Remapper>::iterator r = remappers.begin();
	  r != remappers.end(); ++r ) {
      r->second.finalize();
    }

  }

  theRemappersFilled = true;

}

void BinSampler::runIteration(unsigned long points, bool progress) {

  progress_display* progressBar = 0;
  if ( progress ) {
    Repository::clog() << "integrating " << ANSI::red << process()<< ANSI::reset << ", iteration "
		       << (iterations().size() + 1);
    progressBar = new progress_display(points,Repository::clog());
  }

  double w=0.;
  double maxweight=0;
  int numlastmax=0;
  unsigned long countzero =0;
  int newmax=0;
  for ( unsigned long k = 0; k < points; ++k,++countzero ) {
    if (countzero>=points)break;  
    w=abs(generate());
    
    if(theNonZeroInPresampling && w==0.0){
      k--;
      continue;
    }
    if (w!=0.0)
      countzero =0;
    numlastmax++;
    if (theHalfPoints&&maxweight<w&&
        numlastmax<(int)(points/2.)){
      if(++newmax>theMaxNewMax){
          throw Exception()
             << "\n---------------------------------------------------\n\n"
             << "To many new Maxima.\n\n"
             << "* With the option:\n\n"
             << "* set Sampler:BinSampler:HalfPoints Yes\n\n"
             << "* for every new maximum weight found until the half of the persampling points\n"
             << "* the counter is set to zero. We count the number of new maxima.\n" 
             << "* You have reached: "<<newmax<<"\n"
             << "* Did you apply reasonable cuts to the process?\n"
             << "* You can set the maximum allowed new maxima by:"
             << "* set Sampler:BinSampler:MaxNewMax N\n\n"  
             << "---------------------------------------------------\n"
             << Exception::abortnow;
      }


      maxweight=w;
      k=0;
      numlastmax=0;
    }

    if ( progress ) {
      ++(*progressBar);
    }

  }

  if ( progress ) {
    Repository::clog() << "integrated ( " << ANSI::yellow 
		       << averageWeight() << " +/- " << sqrt(averageWeightVariance())
		       << ANSI::reset << " ) nb\nepsilon = "
		       << (abs(maxWeight()) != 0. ? averageAbsWeight()/abs(maxWeight()) : 0.);
    if ( !iterations().empty() )
      Repository::clog() << " chi2 = " << chi2();
    Repository::clog() << "\n";
    Repository::clog() << "---------------------------------------------------\n";
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
  if ( !integrated() ) {
    unsigned long points = initialPoints();
    for ( unsigned long k = 0; k < nIterations(); ++k ) {
      runIteration(points,progress);
      if ( k < nIterations() - 1 ) {
	points = (unsigned long)(points*enhancementFactor());
	adapt();
	nextIteration();
      }
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
  double sumofweights2=0.;
  for ( map<double,double >::const_iterator b = binsw1.begin();b != binsw1.end(); ++b )
       sumofweights2+=b->second;  
  map<double,double >::const_iterator b2 = binsw1.begin();
  if ( sumofweights == 0 ) {
    cerr << "Not enough statistic accumulated for "
	 << process << " skipping random number diagnostic.\n"
	 << flush;
    return;
  }

  for ( map<double,double >::const_iterator b = bins.begin();
	b != bins.end(); ++b, ++b2) {
      out << " " << b->first
	  << " " << b->second/sumofweights*100.
	  << " " << b2->second/sumofweights2*100.
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
    << ".dat' u ($1):($2)  w boxes  lc rgbcolor \"blue\" t '{\\tiny "<<process <<" }',";
    gpout << " '" << prefix2+fname.str();
    gpout << ".dat' u ($1):($3)  w boxes  lc rgbcolor \"red\" t '';";
  gpout << "reset\n";
}








// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void BinSampler::persistentOutput(PersistentOStream & os) const {
  MultiIterationStatistics::put(os);
  os << theBias << theWeighted << theInitialPoints << theNIterations 
     << theEnhancementFactor << theNonZeroInPresampling << theHalfPoints 
     << theMaxNewMax << theReferenceWeight
     << theBin << theInitialized << theLastPoint
     << theEventHandler << theSampler << theRandomNumbers
     << theRemapperPoints << theRemapChannelDimension
     << theLuminosityMapperBins << theGeneralMapperBins << theKappa;
}

void BinSampler::persistentInput(PersistentIStream & is, int) {
  MultiIterationStatistics::get(is);
  is >> theBias >> theWeighted >> theInitialPoints >> theNIterations 
     >> theEnhancementFactor >> theNonZeroInPresampling >> theHalfPoints 
     >> theMaxNewMax >> theReferenceWeight
     >> theBin >> theInitialized >> theLastPoint
     >> theEventHandler >> theSampler >> theRandomNumbers
     >> theRemapperPoints >> theRemapChannelDimension
     >> theLuminosityMapperBins >> theGeneralMapperBins >> theKappa;
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

  static Switch<BinSampler,bool> interfaceNonZeroInPresampling
    ("NonZeroInPresampling",
     "Switch on to count only non zero weights in presampling.",
     &BinSampler::theNonZeroInPresampling, true, false, false);
  static SwitchOption interfaceNonZeroInPresamplingYes
    (interfaceNonZeroInPresampling,
     "Yes",
     "",
     true);
  static SwitchOption interfaceNonZeroInPresamplingNo
    (interfaceNonZeroInPresampling,
     "No",
     "",
     false);

  static Switch<BinSampler,bool> interfaceHalfPoints
    ("HalfPoints",
     "Switch on to reset the counter of points if new maximumis was found in the first 1/2 points.",
     &BinSampler::theHalfPoints, true, false, false);
  static SwitchOption interfaceHalfPointsYes
    (interfaceHalfPoints,
     "Yes",
     "",
     true);
  static SwitchOption interfaceHalfPointsNo
    (interfaceHalfPoints,
     "No",
     "",
     false);

  static Parameter<BinSampler,int> interfaceMaxNewMax
    ("MaxNewMax",
     "The maximum number of allowed new maxima in combination with the HalfPoints option.",
     &BinSampler::theMaxNewMax, 30, 1, 0,
     false, false, Interface::lowerlim);


  static Parameter<BinSampler,string> interfaceRandomNumbers
    ("RandomNumbers",
     "Prefix for distributions of the random numbers.",
     &BinSampler::theRandomNumbers, "",
     false, false);

  static Parameter<BinSampler,unsigned long> interfaceRemapperPoints
    ("RemapperPoints",
     "The number of points to be used for filling remappers.",
     &BinSampler::theRemapperPoints, 10000, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<BinSampler,bool> interfaceRemapChannelDimension
    ("RemapChannelDimension",
     "Switch on remapping of the channel dimension.",
     &BinSampler::theRemapChannelDimension, true, false, false);
  static SwitchOption interfaceRemapChannelDimensionYes
    (interfaceRemapChannelDimension,
     "Yes",
     "",
     true);
  static SwitchOption interfaceRemapChannelDimensionNo
    (interfaceRemapChannelDimension,
     "No",
     "",
     false);

  static Parameter<BinSampler,unsigned long> interfaceLuminosityMapperBins
    ("LuminosityMapperBins",
     "The number of bins to be used for remapping parton luminosities.",
     &BinSampler::theLuminosityMapperBins, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<BinSampler,unsigned long> interfaceGeneralMapperBins
    ("GeneralMapperBins",
     "The number of bins to be used for remapping other phase space dimensions.",
     &BinSampler::theGeneralMapperBins, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<BinSampler,double> interfaceRemapperMinSelection
    ("RemapperMinSelection",
     "The minimum bin selection probability for remappers.",
     &BinSampler::theRemapperMinSelection, 0.00001, 0.0, 1.0,
     false, false, Interface::limited);


  static Parameter<BinSampler,double> interfaceKappa
    ("Kappa",
     "In AllmostUnweighted mode unweight to Kappa ReferenceWeight.",
     &BinSampler::theKappa, 1., 0.000001, 1.0,
     false, false, Interface::limited);



}

