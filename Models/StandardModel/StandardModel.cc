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

namespace Herwig {
using namespace ThePEG;

StandardModel::~StandardModel() {}

void StandardModel::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex <<_theFFPVertex << _theFFGVertex << _theFFWVertex 
     << _theFFHVertex << _theWWHVertex << _theGGGGVertex << _theWWWWVertex
     << _theGGGVertex << _theWWWVertex << _theRunningMass;
}

void StandardModel::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _theFFGVertex >> _theFFWVertex
     >> _theFFHVertex >> _theWWHVertex >> _theGGGGVertex >> _theWWWWVertex
     >> _theGGGVertex >> _theWWWVertex >> _theRunningMass;
}

ClassDescription<StandardModel> StandardModel::initStandardModel;
// Definition of the static class description member.

void StandardModel::Init() {

static Reference<StandardModel,Herwig::Helicity::FFVVertex> interfaceVertexFFZ
  ("Vertex/FFZ",
   "Reference to the Standard Model FFZ Vertex",
   &StandardModel::_theFFZVertex, false, false, true, false);

static Reference<StandardModel,Herwig::Helicity::FFVVertex> interfaceVertexFFP
  ("Vertex/FFP",
   "Reference to the Standard Model FFP Vertex",
   &StandardModel::_theFFPVertex, false, false, true, false);

static Reference<StandardModel,Herwig::Helicity::FFVVertex> interfaceVertexFFG
  ("Vertex/FFG",
   "Reference to the Standard Model FFG Vertex",
   &StandardModel::_theFFGVertex, false, false, true, false);

static Reference<StandardModel,Herwig::Helicity::FFVVertex> interfaceVertexFFW
  ("Vertex/FFW",
   "Reference to the Standard Model FFW Vertex",
   &StandardModel::_theFFWVertex, false, false, true, false);


static Reference<StandardModel,Herwig::Helicity::FFSVertex> interfaceVertexFFH
  ("Vertex/FFH",
   "Reference to the Standard Model FFH Vertex.",
   &StandardModel::_theFFHVertex, false, false, true, false);

static Reference<StandardModel,Herwig::Helicity::VVVVertex> interfaceVertexGGG
  ("Vertex/GGG",
   "Reference to the Standard Model GGG Vertex",
   &StandardModel::_theGGGVertex, false, false, true, false, false);

static Reference<StandardModel,Herwig::Helicity::VVVVertex> interfaceVertexWWW
  ("Vertex/WWW",
   "Reference to the Standard Model WWW Vertex",
   &StandardModel::_theWWWVertex, false, false, true, false, false);


static Reference<StandardModel,Herwig::Helicity::VVSVertex> interfaceVertexWWH
  ("Vertex/WWH",
   "Reference to the Standard Model WWH Vertex",
   &StandardModel::_theWWHVertex, false, false, true, false);


static Reference<StandardModel,Herwig::Helicity::VVVVVertex> interfaceVertexWWWW
  ("Vertex/WWWW",
   "Reference to the Standard Model WWWW Vertex",
   &StandardModel::_theWWWWVertex, false, false, true, false);

static Reference<StandardModel,Herwig::Helicity::VVVVVertex> interfaceVertexGGGG
  ("Vertex/GGGG",
   "Reference to the Standard Model GGGG Vertex",
   &StandardModel::_theGGGGVertex, false, false, true, false);


static Reference<StandardModel,RunningMassBase> interfaceRunningMass
  ("RunningMass",
   "Reference to the running mass object",
   &StandardModel::_theRunningMass, false, false, true, false);

static ClassDocumentation<StandardModel> documentation
  ("The \\classname{StandardModel} class inherits from StandardModelBase"
   "and supplies additional couplings and access to the StandardModel"
   "vertices for helicity amplitude calculations" );
}

}

