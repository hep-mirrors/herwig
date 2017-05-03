// -*- C++ -*-
//
// MonacoSampler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MonacoSampler class.
//

#include "MonacoSampler.h"
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

#include "MonacoSampler.h"
#include "Herwig/Sampling/GeneralSampler.h"

using namespace Herwig;

MonacoSampler::MonacoSampler() 
  : BinSampler(),
    theAlpha(0.875),
    theGridDivisions(48),
    theIterationPoints(0) {}

MonacoSampler::~MonacoSampler() {}

IBPtr MonacoSampler::clone() const {
  return new_ptr(*this);
}

IBPtr MonacoSampler::fullclone() const {
  return new_ptr(*this);
}

double MonacoSampler::generate() {
  double w = 1.;
//  cout<<"\npoint: ";
  std::valarray<int> upperb(dimension());
  for ( int k = 0; k < dimension(); ++k ) {
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
    lastPoint()[k] = glower + (div-upperb[k])*gdiff;
    w *= gdiff * theGridDivisions;
  }
//    cout<<lastPoint()[k]<<" ";
  
  try {
    w *= eventHandler()->dSigDR(lastPoint()) / nanobarn;
  } catch (Veto&) {
    w = 0.0;
  } catch (...) {
    throw;
  }
// only store numbers
  double wgt = w;
  if ( ! isfinite(wgt) ) wgt = 0;
// save results for later grid optimization 
  theIterationPoints++;
  for ( int k = 0; k < dimension(); ++k ) {
    theGridData(k,upperb[k]) += wgt*wgt;
  }

  if (randomNumberString()!="") 
  for ( size_t k = 0; k < lastPoint().size(); ++k ) {
    RandomNumberHistograms[RandomNumberIndex(id(),k)].first.book(lastPoint()[k],wgt);
    RandomNumberHistograms[RandomNumberIndex(id(),k)].second+=wgt;
  }

  if ( !weighted() && initialized() ) {
    double p = min(abs(w),kappa()*referenceWeight())/(kappa()*referenceWeight());
    double sign = w >= 0. ? 1. : -1.;
    if ( p < 1 && UseRandom::rnd() > p )
      w = 0.;
    else
      w = sign*max(abs(w),kappa()*referenceWeight());
  }
  select(w);
  assert(kappa()==1.||sampler()->almostUnweighted());
  if ( w != 0.0 )
    accept();
  return w;
}

void MonacoSampler::saveGrid() const {
  XML::Element grid = toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);
}

bool MonacoSampler::existsGrid() const {
  list<XML::Element>::iterator git = sampler()->grids().children().begin();
  for ( ; git != sampler()->grids().children().end(); ++git ) {
    if ( git->type() != XML::ElementTypes::Element )
      continue;
    if ( git->name() != "Monaco" )
      continue;
    string proc;
    git->getFromAttribute("process",proc);
    if ( proc == id() ) 
      return true;
  }
  return false;
}

void MonacoSampler::initialize(bool progress) {

//read in grid
  bool haveGrid = false;
  list<XML::Element>::iterator git = sampler()->grids().children().begin();
  for ( ; git != sampler()->grids().children().end(); ++git ) {
    if ( git->type() != XML::ElementTypes::Element )
      continue;
    if ( git->name() != "Monaco" )
      continue;
    string proc;
    git->getFromAttribute("process",proc);
    if ( proc == id() ) {
      haveGrid = true;
      break;
    }
  }

  if ( haveGrid ) {
    fromXML(*git);
    sampler()->grids().erase(git);
    didReadGrids();
  } else {
// flat grid
    theGrid.resize(dimension(),theGridDivisions);
    for (int k = 0; k < dimension(); k++) 
      for (size_t l = 0; l < theGridDivisions; l++)  
        theGrid(k,l) = (l+1)/static_cast<double>(theGridDivisions);
    theGridData = boost::numeric::ublas::zero_matrix<double>(dimension(),theGridDivisions);
    theIterationPoints = 0;
  }

  lastPoint().resize(dimension());
  if (randomNumberString()!="") 
  for(size_t i=0;i<lastPoint().size();i++){
     RandomNumberHistograms[RandomNumberIndex(id(),i)] = make_pair( RandomNumberHistogram(),0.);
  }

  if ( initialized() ) {
    if ( !hasGrids() )
      throw Exception() << "MonacoSampler: Require existing grid when starting to run.\n"
			<< "Did you miss setting --setupfile?"
			<< Exception::abortnow;
    return;
  }

  if ( haveGrid ) {
    if ( !integrated() ) {
      runIteration(initialPoints(),progress);
      adapt();
    }
    isInitialized();
    return;
  }

//   if ( !sampler()->grids().children().empty() ) {
//     nIterations(1);
//   }
  unsigned long points = initialPoints();
  for ( unsigned long k = 0; k < nIterations(); ++k ) {
    runIteration(points,progress);
    if ( k < nIterations() - 1 ) {
      points = (unsigned long)(points*enhancementFactor());
      adapt();
      nextIteration();
    }
  }
  adapt();
  didReadGrids();
  isInitialized();
}

void MonacoSampler::adapt() {

  int dim = dimension();

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

  theGridData = boost::numeric::ublas::zero_matrix<double>(dimension(),theGridDivisions);
  theIterationPoints = 0;
}

void MonacoSampler::finalize(bool) {
// save grid
  adapt();
  XML::Element grid = MonacoSampler::toXML();
  grid.appendAttribute("process",id());
  sampler()->grids().append(grid);

  if (randomNumberString()!="")  
    for ( map<RandomNumberIndex,pair<RandomNumberHistogram,double> >::
	  const_iterator b = RandomNumberHistograms.begin();
	b != RandomNumberHistograms.end(); ++b ) {
       b->second.first.dump(randomNumberString(), b->first.first,shortprocess(),b->first.second);
  }

}

void MonacoSampler::fromXML(const XML::Element& grid) {
  int dim = 0;
  grid.getFromAttribute("Dimension",dim);
  if ( dim != dimension() ) {
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

  theGridData = boost::numeric::ublas::zero_matrix<double>(dimension(),theGridDivisions);

}

XML::Element MonacoSampler::toXML() const {

  XML::Element grid(XML::ElementTypes::Element,"Monaco");
  grid.appendAttribute("Dimension",dimension());
  grid.appendAttribute("GridDivisions",theGridDivisions);

  for ( int k = 0; k < dimension(); ++k ) {
    XML::Element gridvector(XML::ElementTypes::Element,"GridVector");
    gridvector.appendAttribute("Index",k);
    ostringstream bdata;
    bdata << setprecision(17);
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


void MonacoSampler::persistentOutput(PersistentOStream & os) const {
  BinSampler::put(os);
  os << theAlpha << theGridDivisions;
}

void MonacoSampler::persistentInput(PersistentIStream & is, int) {
  BinSampler::get(is);
  is >> theAlpha >> theGridDivisions;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MonacoSampler,BinSampler>
  describeHerwigMonacoSampler("Herwig::MonacoSampler", "HwSampling.so");

void MonacoSampler::Init() {

  static ClassDocumentation<MonacoSampler> documentation
    ("MonacoSampler samples XCombs bins. This implementation performs weighted MC integration using Monaco, an adapted Vegas algorithm.");

  static Parameter<MonacoSampler,double> interfaceAlpha
    ("Alpha",
     "Rate of grid modification (0 for no modification).",
     &MonacoSampler::theAlpha, 0.875, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MonacoSampler,size_t> interfaceGridDivisions
    ("GridDivisions",
     "The number of divisions per grid dimension.",
     &MonacoSampler::theGridDivisions, 48, 1, 0,
     false, false, Interface::lowerlim);

}
