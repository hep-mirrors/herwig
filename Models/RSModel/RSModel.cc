// -*- C++ -*-
//
// RSModel.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModel class.
//

#include "RSModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

RSModel::RSModel() : _theLambda_pi(10000*GeV) {}

void RSModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(_theLambda_pi,GeV) 
     << _theFFGRVertex << _theVVGRVertex << _theSSGRVertex 
     << _theFFVGRVertex << _theVVVGRVertex;
}

void RSModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_theLambda_pi,GeV) 
     >> _theFFGRVertex >> _theVVGRVertex >> _theSSGRVertex
     >> _theFFVGRVertex >> _theVVVGRVertex;
}

ClassDescription<RSModel> RSModel::initRSModel;
// Definition of the static class description member.

void RSModel::Init() {
  

static Reference<RSModel,ThePEG::Helicity::AbstractFFTVertex> interfaceVertexFFGR
  ("Vertex/FFGR",
   "Reference to the fermion-fermion-graviton vertex",
   &RSModel::_theFFGRVertex, false, false, true, false, false);

static Reference<RSModel,ThePEG::Helicity::AbstractVVTVertex> interfaceVertexVVGR
  ("Vertex/VVGR",
   "Reference to the vector-vector-graviton vertex",
   &RSModel::_theVVGRVertex, false, false, true, false, false);

static Reference<RSModel,ThePEG::Helicity::AbstractSSTVertex> interfaceVertexSSGR
  ("Vertex/SSGR",
   "Reference to the scalar-scalar-graviton vertex",
   &RSModel::_theSSGRVertex, false, false, true, false, false);

static Reference<RSModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFVGR
  ("Vertex/FFVGR",
   "Reference to the fermion-antifermion-vector graviton vertex",
   &RSModel::_theFFVGRVertex, false, false, true, false, false);
static Reference<RSModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexVVVGR
  ("Vertex/VVVGR",
   "Reference to the three vector graviton vertex",
   &RSModel::_theVVVGRVertex, false, false, true, false, false);
  
static Parameter<RSModel,Energy> interfaceLambda_pi
  ("Lambda_pi",
   "The coupling of the graviton to matter",
   &RSModel::_theLambda_pi, GeV, 10000*GeV, ZERO, 1.0e12*GeV,
   false, false, false);

  static ClassDocumentation<RSModel> documentation
    ("The RSModel class replaces the Standard Model class for the"
     " RS model");
  
}

}
