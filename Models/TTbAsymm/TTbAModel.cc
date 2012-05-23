
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModel class.
//

#include "TTbAModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;


void TTbAModel::doinit()  {
  addVertex(_theWPTDVertex);
  addVertex(_theZPQQVertex);
  addVertex(_theAGQQVertex);
  addVertex(_theSU2XVertex);
  StandardModel::doinit();
}

TTbAModel::TTbAModel(): _gWPTD_L(1.0), _gWPTD_R(1.0),_gZPTU_L(1.0), _gZPTU_R(1.0),_gZPUU_L(1.0), _gZPUU_R(1.0),_gZPCC_L(1.0), _gZPCC_R(1.0),_gAGQQ_L(1.0), _gAGQQ_R(1.0),_gAGTT_L(1.0), _gAGTT_R(1.0), _alphaXparam(0.060), _costhetaXparam(0.95), _modelselect(1) {}

IBPtr TTbAModel::clone() const {
  return new_ptr(*this);
}
IBPtr TTbAModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TTbAModel::persistentOutput(PersistentOStream & os) const {
  os << _theWPTDVertex
     << _theZPQQVertex
     << _theAGQQVertex
     << _theSU2XVertex
     << _gWPTD_L
     << _gWPTD_R
     << _gZPTU_L
     << _gZPTU_R
     << _gZPUU_L
     << _gZPUU_R
     << _gZPCC_L
     << _gZPCC_R
     << _gAGQQ_L
     << _gAGQQ_R
     << _gAGTT_L
     << _gAGTT_R
     << _modelselect;
}

void TTbAModel::persistentInput(PersistentIStream & is, int) {
  is >> _theWPTDVertex
     >> _theZPQQVertex
     >> _theAGQQVertex
     >> _theSU2XVertex
     >> _gWPTD_L
     >> _gWPTD_R
     >> _gZPTU_L
     >> _gZPTU_R
     >> _gZPUU_L
     >> _gZPUU_R
     >> _gZPCC_L
     >> _gZPCC_R
     >> _gAGQQ_L
     >> _gAGQQ_R
     >> _gAGTT_L
     >> _gAGTT_R
     >> _modelselect;
}

ClassDescription<TTbAModel> TTbAModel::initTTbAModel;
// Definition of the static class description member.

void TTbAModel::Init() {
  
  static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexWPTD
  ("Vertex/WPTD",
   "Reference to the W prime Top Down vertex",
   &TTbAModel::_theWPTDVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexZPQQ
  ("Vertex/ZPQQ",
   "Reference to the Z prime Quark-Antiquark vertex",
   &TTbAModel::_theZPQQVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexAGQQ
  ("Vertex/AGQQ",
   "Reference to the Axigluon Quark-Antiquark vertex",
   &TTbAModel::_theAGQQVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexSU2X
  ("Vertex/SU2X",
   "Reference to the non-Abelian SU(2)_X vertex",
   &TTbAModel::_theSU2XVertex, false, false, true, false, false);




  static Parameter<TTbAModel, double> interfaceWPTDLCoupling
    ("WPTDLCoupling",
     "The left-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceWPTDRCoupling
    ("WPTDRCoupling",
     "The right-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceAGQQLCoupling
    ("AGQQLCoupling",
     "The left-handed axigluon coupling to q-qbar",
     &TTbAModel::_gAGQQ_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

static Parameter<TTbAModel, double> interfaceAGQQRCoupling
    ("AGQQRCoupling",
     "The right-handed axigluon coupling to q-qbar",
     &TTbAModel::_gAGQQ_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);


  static Parameter<TTbAModel, double> interfaceAGTTLCoupling
    ("AGTTLCoupling",
     "The left-handed axigluon coupling to t-tbar",
     &TTbAModel::_gAGTT_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

static Parameter<TTbAModel, double> interfaceAGTTRCoupling
    ("AGTTRCoupling",
     "The right-handed axigluon coupling to t-tbar",
     &TTbAModel::_gAGTT_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPTULCoupling
    ("ZPTULCoupling",
     "The left-handed Z prime coupling to top up",
     &TTbAModel::_gZPTU_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPTURCoupling
    ("ZPTURCoupling",
     "The right-handed Z prime coupling to top up",
     &TTbAModel::_gZPTU_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPUULCoupling
    ("ZPUULCoupling",
     "The left-handed Z prime coupling to up upbar",
     &TTbAModel::_gZPUU_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPUURCoupling
    ("ZPUURCoupling",
     "The right-handed Z prime coupling to up upbar",
     &TTbAModel::_gZPUU_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPCCLCoupling
    ("ZPCCLCoupling",
     "The left-handed Z prime coupling to char charmbar",
     &TTbAModel::_gZPCC_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPCCRCoupling
    ("ZPCCRCoupling",
     "The right-handed Z prime coupling to charm charmbar",
     &TTbAModel::_gZPCC_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);


  static Parameter<TTbAModel, double> interfaceSU2Xcostheta
    ("SU2Xcostheta",
     "Misalignment parameter of SU(2)_X model",
     &TTbAModel::_costhetaXparam, 0.95, -1.0, 1.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceSU2Xalpha
    ("SU2Xalpha",
     "alphaX coupling constant",
     &TTbAModel::_alphaXparam, 0.060, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, int> interfacemodelselect
    ("modelselect",
     "Selet which model to run",
     &TTbAModel::_modelselect, 0, 0, 4,
     false, false, Interface::limited);


  static ClassDocumentation<TTbAModel> documentation
    ("There is no documentation for the TTbAModel class");

}

