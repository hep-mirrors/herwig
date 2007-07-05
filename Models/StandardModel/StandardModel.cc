// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModel class.
//

#include "StandardModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/General/ModelGenerator.h"

namespace Herwig {
using namespace ThePEG;

IBPtr StandardModel::clone() const {
  return new_ptr(*this);
}

IBPtr StandardModel::fullclone() const {
  return new_ptr(*this);
}

StandardModel::StandardModel() {}

StandardModel::~StandardModel() {}

StandardModel::StandardModel(const StandardModel & x)
  : StandardModelBase(x) 
  ,_theFFZVertex(x._theFFZVertex)
  ,_theFFPVertex(x._theFFPVertex) ,_theFFGVertex(x._theFFGVertex)
  ,_theFFWVertex(x._theFFWVertex) ,_theFFHVertex(x._theFFHVertex)
  ,_theWWHVertex(x._theWWHVertex) 
  ,_theGGGVertex(x._theGGGVertex) 
  ,_theWWWVertex(x._theWWWVertex) , _theGGGGVertex(x._theGGGGVertex) 
  ,_theWWWWVertex(x._theWWWWVertex),_theHGGVertex(x._theHGGVertex)
 ,_theHPPVertex(x._theHPPVertex) ,_vertexlist(x._vertexlist)
  ,_theRunningMass(x._theRunningMass),_theModelGenerator(x._theModelGenerator) {}

void StandardModel::doinit() throw(InitException) {
  StandardModelBase::doinit();
  if(_theRunningMass) {
    _theRunningMass->init();
  }
  //add Standard Model vertices
  addVertex(_theFFZVertex);
  addVertex(_theFFPVertex);
  addVertex(_theFFGVertex);
  addVertex(_theFFWVertex);
  addVertex(_theFFHVertex);
  addVertex(_theWWHVertex);
  addVertex(_theGGGVertex);
  addVertex(_theWWWVertex);
  addVertex(_theGGGGVertex);
  addVertex(_theWWWWVertex);
  addVertex(_theHGGVertex);
  addVertex(_theHPPVertex);
}

void StandardModel::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex <<_theFFPVertex << _theFFGVertex << _theFFWVertex 
     << _theFFHVertex << _theWWHVertex << _theGGGGVertex << _theWWWWVertex
     << _theGGGVertex << _theWWWVertex  << _theHGGVertex  << _theHPPVertex 
     << _theRunningMass << _vertexlist << _theModelGenerator;
}

void StandardModel::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _theFFGVertex >> _theFFWVertex
     >> _theFFHVertex >> _theWWHVertex >> _theGGGGVertex >> _theWWWWVertex
     >> _theGGGVertex >> _theWWWVertex >> _theHGGVertex  >> _theHPPVertex 
     >> _theRunningMass >> _vertexlist >> _theModelGenerator;
}

ClassDescription<StandardModel> StandardModel::initStandardModel;
// Definition of the static class description member.

void StandardModel::Init() {

static Reference<StandardModel,ThePEG::Helicity::FFVVertex> interfaceVertexFFZ
  ("Vertex/FFZ",
   "Reference to the Standard Model FFZ Vertex",
   &StandardModel::_theFFZVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::FFVVertex> interfaceVertexFFP
  ("Vertex/FFP",
   "Reference to the Standard Model FFP Vertex",
   &StandardModel::_theFFPVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::FFVVertex> interfaceVertexFFG
  ("Vertex/FFG",
   "Reference to the Standard Model FFG Vertex",
   &StandardModel::_theFFGVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::FFVVertex> interfaceVertexFFW
  ("Vertex/FFW",
   "Reference to the Standard Model FFW Vertex",
   &StandardModel::_theFFWVertex, false, false, true, false);


static Reference<StandardModel,ThePEG::Helicity::FFSVertex> interfaceVertexFFH
  ("Vertex/FFH",
   "Reference to the Standard Model FFH Vertex.",
   &StandardModel::_theFFHVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::VVVVertex> interfaceVertexGGG
  ("Vertex/GGG",
   "Reference to the Standard Model GGG Vertex",
   &StandardModel::_theGGGVertex, false, false, true, false, false);

static Reference<StandardModel,ThePEG::Helicity::VVVVertex> interfaceVertexWWW
  ("Vertex/WWW",
   "Reference to the Standard Model WWW Vertex",
   &StandardModel::_theWWWVertex, false, false, true, false, false);


static Reference<StandardModel,ThePEG::Helicity::VVSVertex> interfaceVertexWWH
  ("Vertex/WWH",
   "Reference to the Standard Model WWH Vertex",
   &StandardModel::_theWWHVertex, false, false, true, false);


static Reference<StandardModel,ThePEG::Helicity::VVVVVertex> interfaceVertexWWWW
  ("Vertex/WWWW",
   "Reference to the Standard Model WWWW Vertex",
   &StandardModel::_theWWWWVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::VVVVVertex> interfaceVertexGGGG
  ("Vertex/GGGG",
   "Reference to the Standard Model GGGG Vertex",
   &StandardModel::_theGGGGVertex, false, false, true, false);

static Reference<StandardModel,ThePEG::Helicity::GeneralSVVVertex> interfaceVertexHGG
  ("Vertex/HGG",
   "Reference to the StandardModel HGG Vertex",
   &StandardModel::_theHGGVertex, false, false, true, false);
 
static Reference<StandardModel,ThePEG::Helicity::GeneralSVVVertex> interfaceVertexHPP
  ("Vertex/HPP",
   "Reference to StandardModel HPPVertex",
   &StandardModel::_theHPPVertex, false, false, true, false);

static Reference<StandardModel,RunningMassBase> interfaceRunningMass
  ("RunningMass",
   "Reference to the running mass object",
   &StandardModel::_theRunningMass, false, false, true, false);

static Reference<StandardModel,Herwig::ModelGenerator> interfaceModelGenerator
     ("ModelGenerator",
      "Pointer to ModelGenerator class",
      &StandardModel::_theModelGenerator, false, false, true, true);

static ClassDocumentation<StandardModel> documentation
  ("The StandardModel class inherits from StandardModelBase"
   "and supplies additional couplings and access to the StandardModel"
   "vertices for helicity amplitude calculations" );

}
}

