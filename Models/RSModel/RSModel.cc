// -*- C++ -*-
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

RSModel::~RSModel() {}

void RSModel::persistentOutput(PersistentOStream & os) const {
  os << _theLambda_pi << _theFFGRVertex << _theVVGRVertex << _theSSGRVertex 
     << _theFFVGRVertex << _theVVVGRVertex;
}

void RSModel::persistentInput(PersistentIStream & is, int) {
  is >> _theLambda_pi >> _theFFGRVertex >> _theVVGRVertex >> _theSSGRVertex
     >> _theFFVGRVertex >> _theVVVGRVertex;
}

ClassDescription<RSModel> RSModel::initRSModel;
// Definition of the static class description member.

void RSModel::Init() {
  

static Reference<RSModel,Herwig::Helicity::FFTVertex> interfaceVertexFFGR
  ("Vertex/FFGR",
   "Reference to the fermion-fermion-graviton vertex",
   &RSModel::_theFFGRVertex, false, false, true, false, false);

static Reference<RSModel,Herwig::Helicity::VVTVertex> interfaceVertexVVGR
  ("Vertex/VVGR",
   "Reference to the vector-vector-graviton vertex",
   &RSModel::_theVVGRVertex, false, false, true, false, false);

static Reference<RSModel,Herwig::Helicity::SSTVertex> interfaceVertexSSGR
  ("Vertex/SSGR",
   "Reference to the scalar-scalar-graviton vertex",
   &RSModel::_theSSGRVertex, false, false, true, false, false);

static Reference<RSModel,Herwig::Helicity::FFVTVertex> interfaceVertexFFVGR
  ("Vertex/FFVGR",
   "Reference to the fermion-antifermion-vector graviton vertex",
   &RSModel::_theFFVGRVertex, false, false, true, false, false);

static Reference<RSModel,Herwig::Helicity::VVVTVertex> interfaceVertexVVVGR
  ("Vertex/VVVGR",
   "Reference to the three vector graviton vertex",
   &RSModel::_theVVVGRVertex, false, false, true, false, false);
  
static Parameter<RSModel,Energy> interfaceLambda_pi
  ("Lambda_pi",
   "The coupling of the graviton to matter",
   &RSModel::_theLambda_pi, GeV, 10000*GeV, 0*GeV, 1.0e12*GeV,
   false, false, false);

  static ClassDocumentation<RSModel> documentation
    ("The \\classname{RSModel} class replaces the Standard Model class for the"
     " RS model");
  
}

}
