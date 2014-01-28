// -*- C++ -*-
//
// combineRuns.cpp is a part of myStatistics
// Copyright (C) 2012-2013 Simon Platzer
//
// myStatistics is licenced under version 2 of the GPL, see COPYING for details.
//

#include "Run.h"

#include "Herwig++/Utilities/XML/ElementIO.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace Statistics;
using namespace XML;
using namespace std;

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
  combined.fromXML(ElementIO::get(firstIn));
  combined.name(outName);

  for ( int k = 3; k < argc; ++k ) {
    Run next;
    ifstream nextIn(argv[k]);
    next.fromXML(ElementIO::get(nextIn));
    combined += next;
  }

  ofstream out(outFileName.c_str());
  ElementIO::put(combined.toXML(),out);

}
