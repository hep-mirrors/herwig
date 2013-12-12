// -*- C++ -*-
//
// CellGridSampler.cpp is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/Utilities/XML/ElementIO.h"

using namespace Herwig;

CellGridSampler::CellGridSampler() {}

CellGridSampler::~CellGridSampler() {}

IBPtr CellGridSampler::clone() const {
  return new_ptr(*this);
}

IBPtr CellGridSampler::fullclone() const {
  return new_ptr(*this);
}

void CellGridSampler::initialize() {
  GeneralSampler::initialize();
  writeGrids();
}

void CellGridSampler::doinit() {
  readGrids();
  GeneralSampler::doinit();
}

void CellGridSampler::dofinish() {
  GeneralSampler::dofinish();
  writeGrids();
}

void CellGridSampler::doinitrun() {
  readGrids();
  GeneralSampler::doinitrun();
}

void CellGridSampler::writeGrids() const {
  string dataName = generator()->filename() + "-grids.xml";
  ofstream out(dataName.c_str());
  XML::ElementIO::put(theGrids,out);
}

void CellGridSampler::readGrids() {
  string dataName = generator()->filename() + "-grids.xml";
  ifstream in(dataName.c_str());
  if ( !in ) {
    theGrids = XML::Element(XML::ElementTypes::Element,"Grids");
    return;
  }
  theGrids = XML::ElementIO::get(in);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void CellGridSampler::persistentOutput(PersistentOStream &) const {}

void CellGridSampler::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<CellGridSampler,GeneralSampler>
  describeHerwigCellGridSampler("Herwig::CellGridSampler", "HwExsample2.so");

void CellGridSampler::Init() {

  static ClassDocumentation<CellGridSampler> documentation
    ("There is no documentation for the CellGridSampler class");

}

