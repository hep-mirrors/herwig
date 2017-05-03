// -*- C++ -*-
//
// makeDistributions.cpp is a part of myStatistics
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//

#include "Run.h"
#include "Distribution.h"

#include "Herwig/Utilities/XML/ElementIO.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Statistics;
using namespace XML;
using namespace std;

template<class T>
inline T sqr(const T& x) {
  return x*x;
}

int main(int argc, char* argv[]) {

  if ( argc < 3 ) {
    cerr << "usage: "
	 << argv[0]
	 << " [outname] [infile 1] { [infile 2] ... [infile n] }\n";
    return 1;
  }

  string outName(argv[1]);
  string outFileName(outName + ".xml");

  Run combined;
  ifstream firstIn(argv[2]);
  firstIn >> setprecision(16);
  XML::Element elemx = ElementIO::get(firstIn);
  combined.fromXML(elemx);
  combined.name(outName);

  for ( int k = 3; k < argc; ++k ) {
    Run next;
    ifstream nextIn(argv[k]);
    nextIn >> setprecision(16);
    elemx = ElementIO::get(nextIn);
    next.fromXML(elemx);
    combined += next;
  }

  double sumOfWeights = combined.sumOfWeights();
  double sumOfSquaredWeights = combined.sumOfSquaredWeights();
  double nPoints = combined.attemptedPoints();

  double integral = sumOfWeights/nPoints;
  double varianceOfIntegral =
    abs(sumOfSquaredWeights/nPoints - sqr(sumOfWeights/nPoints))/(nPoints-1);

  Element elem(ElementTypes::Element,"CrossSections");

  elem.appendAttribute("name",outName);
  elem.appendAttribute("integral",integral);
  elem.appendAttribute("varianceOfIntegral",varianceOfIntegral);

  Element xhistos(XML::ElementTypes::Element,"Distributions");

  for ( map<string,Histogram>::const_iterator h =
	  combined.histograms().begin(); h != combined.histograms().end(); ++h ) {
    Distribution dist(h->second,nPoints);
    // testing only, dump makePlots
    dist.toMakePlots(outName);
    xhistos.append(dist.toXML());
  }

  elem.append(xhistos);

  ofstream out(outFileName.c_str());
  out << setprecision(16);
  ElementIO::put(elem,out);

}
