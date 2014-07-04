// -*- C++ -*-
//
// CellGridSampler2.cpp is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CellGridSampler2 class.
//

#include "CellGridSampler2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <boost/progress.hpp>

#include "CellGridSampler2.h"
#include "Herwig++/Sampling/GeneralSampler.h"

using namespace Herwig;
using namespace ExSample;

CellGridSampler2::CellGridSampler2() 
  : BinSampler(), SimpleCellGrid(),
    theExplorationPoints(1000), theExplorationSteps(8),
    theGain(0.3), theEpsilon(0.01),
    theMinimumSelection(0.0001), theLuminositySplits(0),
    theChannelSplits(0),
    theAlpha(0.875),
    theGridDivisions(48),
    theIterationPoints(0) {}

CellGridSampler2::~CellGridSampler2() {}

IBPtr CellGridSampler2::clone() const {
  return new_ptr(*this);
}

IBPtr CellGridSampler2::fullclone() const {
  return new_ptr(*this);
}


double CellGridSampler2::evaluate(const vector<double>& cellgridlastpoint) {
  double wmonaco = 1.0;
  double w = 1.0;

  std::valarray<int> upperb(monacoDimensions());
  for ( int k = 0; k < monacoDimensions(); ++k ) {
    double div = (1 - UseRandom::rnd()) * theGridDivisions;
    upperb[k] = static_cast<int>(div);
    double gupper, glower;
    if ( upperb[k] <= 0 ) {
      upperb[k] = 0;
      glower = 0.;
      gupper = theGrid(k,0);
    } else if (upperb[k] >= static_cast<int>(theGridDivisions)) {
      upperb[k] = theGridDivisions-1;
      glower = theGrid(k,theGridDivisions-2);
      gupper = theGrid(k,theGridDivisions-1);
    } else {
      glower = theGrid(k,upperb[k]-1);
      gupper = theGrid(k,upperb[k]);
    }
    double gdiff = gupper - glower;
    lastpoint_monaco[k] = glower + (div-upperb[k])*gdiff;
    wmonaco *= gdiff * theGridDivisions;
  }

    int md=0;
    int cgd=0;
    for (int k=0;k < monacoDimensions()+dimension();k++){
      if (theMonacoDimensions.find(k)!=theMonacoDimensions.end()){
	lastPoint()[k]=lastpoint_monaco[md];md++;
      }else{
	lastPoint()[k]=cellgridlastpoint[cgd];cgd++;
      }  
    }
    
    
    
    try {
      w = eventHandler()->dSigDR(lastPoint()) / nanobarn;
    } catch (Veto&) {
      w = 0.0;
    } catch (...) {
      throw;
    }
    
  theIterationPoints++;
  for ( int k = 0; k < monacoDimensions(); ++k ) {
    theGridData(k,upperb[k]) += sqr(w*wmonaco);
  }
    
    
    return w*wmonaco;
  }

double CellGridSampler2::generate() {
  UseRandom rnd;
  
  lastpoint_cellgrid.resize(lastPoint().size()-monacoDimensions());


  pair<double,double> weights = SimpleCellGrid::generate(rnd,*this,lastpoint_cellgrid);

  double w = SimpleCellGrid::integral()*weights.first/weights.second;
  if (randomNumberString()!="") 
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    RandomNumberHistograms[RandomNumberIndex(id(),k)].first.book(lastPoint()[k],abs(w));
    RandomNumberHistograms[RandomNumberIndex(id(),k)].second+=abs(w);

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

void CellGridSampler2::adapt() {
  UseRandom rnd;
  set<SimpleCellGrid*> newCells;
  Monaco_adapt();
  SimpleCellGrid::adapt(theGain,theEpsilon,newCells);
  SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells);
  SimpleCellGrid::setWeights();
  SimpleCellGrid::updateIntegral();
  SimpleCellGrid::minimumSelection(theMinimumSelection);
}

void CellGridSampler2::initialize(bool progress) {
   size_t tmp= theMonacoDimensions.size();
   map<int,double>tmpMD;
   for(map<int,double>::iterator it=theMonacoDimensions.begin(); it!=theMonacoDimensions.end();it++){
     if (it->first<0){
       tmpMD.insert(make_pair(Herwig::BinSampler::dimension()+it->first,1.));      
     }else{
       assert(tmpMD.find(it->first)==tmpMD.end());
       tmpMD.insert(make_pair(it->first,1.));  
     }
   }

   theMonacoDimensions=tmpMD;
   
   assert(tmp== theMonacoDimensions.size());
   
  bool haveCellGridGrid = false;
  list<XML::Element>::iterator gitCG = sampler()->grids().children().begin();
  for ( ; gitCG != sampler()->grids().children().end(); ++gitCG ) {
    if ( gitCG->type() != XML::ElementTypes::Element )
      continue;
    if ( gitCG->name() != "CellGrid" )
      continue;
    string proc;
    gitCG->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveCellGridGrid = true;
      break;
    }
  }
  bool haveMonacoGrid = false;
  list<XML::Element>::iterator gitM = sampler()->grids().children().begin();
  for ( ; gitM != sampler()->grids().children().end(); ++gitM ) {
    if ( gitM->type() != XML::ElementTypes::Element )
      continue;
    if ( gitM->name() != "Monaco" )
      continue;
    string proc;
    gitM->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveMonacoGrid = true;
      break;
    }
  }

  if ( haveCellGridGrid ) {
    SimpleCellGrid::fromXML(*gitCG);
    sampler()->grids().erase(gitCG);
  }
  
  if ( haveMonacoGrid ) {
    Monaco_fromXML(*gitM);
    sampler()->grids().erase(gitM);
  } 
  
  
  
  lastpoint_cellgrid.resize(dimension());
  lastpoint_monaco.resize(monacoDimensions());
  lastPoint().resize(dimension()+monacoDimensions());
  if (randomNumberString()!="") 
  for(size_t i=0;i<lastPoint().size();i++){
     RandomNumberHistograms[RandomNumberIndex(id(),i)] = make_pair( RandomNumberHistogram(),0.);
  }
  
  if ( initialized() ) {
    if ( !haveCellGridGrid &&dimension()!=0)
      throw Exception() << "CellGridSampler2: Require existing CellGrid when starting to run."
			<< Exception::abortnow;
			
    if ( !haveMonacoGrid &&monacoDimensions()!=0)
      throw Exception() << "CellGridSampler2: Require existing MonacoGrid when starting to run."
			<< Exception::abortnow;
    return;
  }

  if ( (haveCellGridGrid||dimension()==0)&&(haveMonacoGrid||monacoDimensions()==0) ) {
    runIteration(initialPoints(),progress);
    isInitialized();
    XML::Element grid = SimpleCellGrid::toXML();
    grid.appendAttribute("process",id());
    sampler()->grids().append(grid);
    
    Monaco_adapt();
    XML::Element Monacogrid = Monaco_toXML();
    Monacogrid.appendAttribute("process",id());
    sampler()->grids().append(Monacogrid);
    return;
  }

  SimpleCellGrid::boundaries(vector<double>(dimension(),0.0),vector<double>(dimension(),1.0));
  SimpleCellGrid::weightInformation().resize(dimension());

//Monacogrid
    theGrid.resize(monacoDimensions(),theGridDivisions);
    for (int k = 0; k < monacoDimensions(); k++) 
      for (size_t l = 0; l < theGridDivisions; l++)  
        theGrid(k,l) = (l+1)/static_cast<double>(theGridDivisions);
    theGridData = boost::numeric::ublas::zero_matrix<double>(monacoDimensions(),theGridDivisions);
    theIterationPoints = 0;
  
  
  
  UseRandom rnd;

  boost::progress_display* progressBar = 0;
  if ( progress ) {
    Repository::clog() << "exploring " << process();
    progressBar = new boost::progress_display(theExplorationSteps,cout);
  }
  std::set<SimpleCellGrid*> newCells;
  
//   if ( pre_adaption_splits().empty() &&
//        (theLuminositySplits || theChannelSplits) ) {
//     const StandardEventHandler& eh = *eventHandler();
//     const StandardXComb& xc = *eh.xCombs()[bin()];
//     the_pre_adaption_splits.resize(dimension(),0);
//     const pair<int,int>& pdims = xc.partonDimensions();
//     if ( theLuminositySplits && dimension() >= pdims.first + pdims.second ) {
//       for ( int n = 0; n < pdims.first; ++n )
// 	the_pre_adaption_splits[n] = theLuminositySplits;
//       for ( int n = dimension() - pdims.second; n < dimension(); ++n )
// 	the_pre_adaption_splits[n] = theLuminositySplits;
//     }
//     if ( theChannelSplits && dimension() > pdims.first + pdims.second ) {
//       the_pre_adaption_splits[pdims.first] = theChannelSplits;
//     }
//   }
//   
//   for(int splitdim=0; splitdim<min(dimension(),(int)pre_adaption_splits().size());splitdim++)
//       SimpleCellGrid::splitter(splitdim,pre_adaption_splits()[splitdim]);
//   
  SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells);

  
  bool notAll = false;
  for ( std::size_t step = 1; step < theExplorationSteps; ++step ) {
    newCells.clear();
    SimpleCellGrid::adapt(theGain,theEpsilon,newCells);
    Monaco_adapt();
    if ( progressBar )
      ++(*progressBar);
    if ( newCells.empty() ) {
      notAll = true;
      break;
    }
    SimpleCellGrid::explore(theExplorationPoints,rnd,*this,newCells);
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
  isInitialized();

  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);
  XML::Element Monacogrid = Monaco_toXML();
  Monacogrid.appendAttribute("process",id());
  sampler()->grids().append(Monacogrid);

}

void CellGridSampler2::finalize(bool) {
  XML::Element grid = SimpleCellGrid::toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);

  Monaco_adapt();
  XML::Element Monacogrid = Monaco_toXML();
  Monacogrid.appendAttribute("process",id());
  sampler()->grids().append(Monacogrid);

  if (randomNumberString()!="")  
    for ( map<RandomNumberIndex,pair<RandomNumberHistogram,double> >::
	  const_iterator b = RandomNumberHistograms.begin();
	b != RandomNumberHistograms.end(); ++b ) {
       b->second.first.dump(randomNumberString(), b->first.first,shortprocess(),b->first.second);
  }

}


void CellGridSampler2::Monaco_adapt() {

  int dim = monacoDimensions();

// refine grid
  std::valarray<double> gridcumul(dim);
  for (int k=0; k<dim; ++k) {
    double gridold = theGridData(k,0);
    double gridnew = theGridData(k,1);
    theGridData(k,0) = (gridold + gridnew) / 2.0;
    gridcumul[k] = theGridData(k,0);
    for (size_t l=1; l<theGridDivisions-1; ++l) {
      theGridData(k,l) = gridold + gridnew;
      gridold = gridnew;
      gridnew = theGridData(k,l+1);
      theGridData(k,l) = (theGridData(k,l) + gridnew) / 3.0;
      gridcumul[k] += theGridData(k,l);
    }
    theGridData(k,theGridDivisions-1) = (gridnew + gridold) / 2.0;
    gridcumul[k] += theGridData(k,theGridDivisions-1);
  }

  for (int k=0; k<dim; ++k) {
     double rc = 0.;
     std::valarray<double> ri(theGridDivisions);
     for (size_t l=0; l<theGridDivisions; ++l) {
        ri[l] = 0.;
        if ((theGridData(k,l) >= 0) && (gridcumul[k] != 0)) {
           theGridData(k,l) = max( 1.0e-30, theGridData(k,l) );
           double gpart = gridcumul[k] / theGridData(k,l);
           ri[l] = pow( (gpart - 1.0) / (gpart * log( gpart )), theAlpha);
        } else {
           ri[l] = pow(  1. / log( 1e30 ), theAlpha);
        }
        rc += ri[l];
     }
     rc /= theGridDivisions;

     double gridold = 0, gridnew = 0.;
     double deltar = 0.;
     unsigned int m = 0;
     std::valarray<double> theGridRowNew(theGridDivisions);

     for (size_t l = 0; l < theGridDivisions; ++l) {
       deltar += ri[l];
       gridold = gridnew;
       gridnew = theGrid(k,l);
       for (; deltar > rc; m++) {
         deltar -= rc;
         theGridRowNew[m] = gridnew - (gridnew - gridold) * deltar / ri[l];
       }
     }
     for (size_t l = 0; l < theGridDivisions-1; ++l) {
        theGrid(k,l) = theGridRowNew[l];
     }
     theGrid(k,theGridDivisions-1) = 1.0;
  }

  theGridData = boost::numeric::ublas::zero_matrix<double>(monacoDimensions(),theGridDivisions);
  theIterationPoints = 0;
}



void CellGridSampler2::Monaco_fromXML(const XML::Element& grid) {
  int dim = 0;
  grid.getFromAttribute("Dimension",dim);
  if ( dim != monacoDimensions() ) {
    throw std::runtime_error("[MonacoSampler] Number of dimensions in grid file does not match expectation.");
  }

  size_t griddivisions = 0;
  grid.getFromAttribute("GridDivisions",griddivisions);
  boost::numeric::ublas::matrix<double> tmpgrid(dim,griddivisions);

  pair<multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator,multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator> cit;
  cit = grid.findAll(XML::ElementTypes::Element,"GridVector");

  if ( cit.first->second == grid.children().end() )
    throw std::runtime_error("[MonacoSampler] Expected a GridVector element.");

  for (multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator iit=cit.first; iit!=cit.second; ++iit) {
    const XML::Element& gridvector = *iit->second;
    

    int k = 0;
    gridvector.getFromAttribute("Index",k);
    if ( k >= dim ) {
      throw std::runtime_error("[MonacoSampler] Index of grid dimension larger than grid size.");
    } else {
      list<XML::Element>::const_iterator git;
      git = gridvector.findFirst(XML::ElementTypes::ParsedCharacterData,"");
      if ( git == gridvector.children().end() )
        throw std::runtime_error("[MonacoSampler] Expected grid data.");
      istringstream bdata(git->content());

      for ( size_t l = 0; l < griddivisions; ++l ) {
        bdata >> tmpgrid(k,l);
      }
    }

  }

// store back into main variable
// if griddivisions do not match, rebin preserving bin density

  theGrid.resize(dim,theGridDivisions);
  theIterationPoints = 0;
  double divratio = griddivisions / static_cast<double>(theGridDivisions);
  for (int k = 0; k < dim; k++) {
    double xold = 0, xnew = 0, deltar = 0;
    size_t l = 0;
    for (size_t m = 0; m < griddivisions; m++) {
      deltar += 1;
      xold = xnew;
      xnew = tmpgrid(k,m);
      for (; deltar > divratio; l++) {
        deltar -= divratio;
        theGrid(k,l) = xnew - (xnew - xold) * deltar;
      }
    }
    theGrid(k,theGridDivisions-1) = 1.0;
  }

  theGridData = boost::numeric::ublas::zero_matrix<double>(monacoDimensions(),theGridDivisions);

}


XML::Element CellGridSampler2::Monaco_toXML() const {

  XML::Element grid(XML::ElementTypes::Element,"Monaco");
  grid.appendAttribute("Dimension",monacoDimensions());
  grid.appendAttribute("GridDivisions",theGridDivisions);

  for ( int k = 0; k < monacoDimensions(); ++k ) {
    XML::Element gridvector(XML::ElementTypes::Element,"GridVector");
    gridvector.appendAttribute("Index",k);
    ostringstream bdata;
    bdata << setprecision(20);
    for ( size_t l = 0; l < theGridDivisions; ++l )
      bdata << theGrid(k,l) << " ";

    XML::Element belem(XML::ElementTypes::ParsedCharacterData,bdata.str());
    gridvector.append(belem);
    grid.append(gridvector);
  }

  return grid;

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void CellGridSampler2::persistentOutput(PersistentOStream & os) const {
  os << theExplorationPoints << theExplorationSteps
     << theGain << theEpsilon << theMinimumSelection
     << the_pre_adaption_splits
     << theLuminositySplits << theChannelSplits << theAlpha << theGridDivisions<<theMonacoDimensions;
}

void CellGridSampler2::persistentInput(PersistentIStream & is, int) {
  is >> theExplorationPoints >> theExplorationSteps
     >> theGain >> theEpsilon >> theMinimumSelection
     >> the_pre_adaption_splits
     >> theLuminositySplits >> theChannelSplits>> theAlpha >> theGridDivisions>>theMonacoDimensions;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<CellGridSampler2,BinSampler>
  describeHerwigCellGridSampler2("Herwig::CellGridSampler2", "HwSampling.so");

void CellGridSampler2::Init() {

  static ClassDocumentation<CellGridSampler2> documentation
    ("CellGridSampler2 samples XCombs bins using CellGrids.");

  static Parameter<CellGridSampler2,size_t> interfaceExplorationPoints
    ("ExplorationPoints",
     "The number of points to use for cell exploration.",
     &CellGridSampler2::theExplorationPoints, 1000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler2,size_t> interfaceExplorationSteps
    ("ExplorationSteps",
     "The number of exploration steps to perform.",
     &CellGridSampler2::theExplorationSteps, 8, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler2,double> interfaceGain
    ("Gain",
     "The gain factor used for adaption.",
     &CellGridSampler2::theGain, 0.3, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<CellGridSampler2,double> interfaceEpsilon
    ("Epsilon",
     "The efficieny threshold used for adaption.",
     &CellGridSampler2::theEpsilon, 0.01, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<CellGridSampler2,double> interfaceMinimumSelection
    ("MinimumSelection",
     "The minimum cell selection probability.",
     &CellGridSampler2::theMinimumSelection, 0.0001, 0.0, 1.0,
     false, false, Interface::limited);
    
  static ParVector<CellGridSampler2,int> interfacethe_pre_adaption_splits
    ("preadaptionsplit",
     "The splittings for each dimension befor adaption.",
     &CellGridSampler2::the_pre_adaption_splits, 1., -1, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler2,int> interfaceLuminositySplits
    ("LuminositySplits",
     "",
     &CellGridSampler2::theLuminositySplits, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler2,int> interfaceChannelSplits
    ("ChannelSplits",
     "",
     &CellGridSampler2::theChannelSplits, 0, 0, 0,
     false, false, Interface::lowerlim);
    
    
  static Parameter<CellGridSampler2,double> interfaceAlpha
    ("Alpha",
     "Rate of grid modification (0 for no modification).",
     &CellGridSampler2::theAlpha, 0.875, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<CellGridSampler2,size_t> interfaceGridDivisions
    ("GridDivisions",
     "The number of divisions per grid dimension.",
     &CellGridSampler2::theGridDivisions, 48, 1, 0,
     false, false, Interface::lowerlim);
    
   static Command<CellGridSampler2> interfacedoMonacoDimension
    ("doMonacoDimension",
     "doMonacoDimension n.",
     &CellGridSampler2::doMonacoDimensions, false);

}

