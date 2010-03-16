// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModel class.
//

#include "LeptoquarkModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void LeptoquarkModel::doinit() {
  addVertex(_theSLQSLQGVertex);
  addVertex(_theSLQSLQGGVertex);
  StandardModel::doinit();
}

LeptoquarkModel::LeptoquarkModel()  {}

LeptoquarkModel::~LeptoquarkModel() {}

IBPtr LeptoquarkModel::clone() const {
  return new_ptr(*this);
}
IBPtr LeptoquarkModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LeptoquarkModel::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os <<  _theSLQSLQGGVertex
     << _theSLQSLQGVertex;

}

void LeptoquarkModel::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _theSLQSLQGGVertex
     >> _theSLQSLQGVertex;
}

ClassDescription<LeptoquarkModel> LeptoquarkModel::initLeptoquarkModel;
// Definition of the static class description member.

void LeptoquarkModel::Init() {
  
  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVSSVertex> interfaceVertexSLQSLQG
  ("Vertex/SLQSLQG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGVertex, false, false, true, false, false);

  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVVSSVertex> interfaceVertexSLQSLQGG
  ("Vertex/SLQSLQGG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGGVertex, false, false, true, false, false);




  static ClassDocumentation<LeptoquarkModel> documentation
    ("There is no documentation for the LeptoquarkModel class");

}

