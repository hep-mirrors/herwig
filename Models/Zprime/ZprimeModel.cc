
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZprimeModel class.
//

#include "ZprimeModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;


void ZprimeModel::doinit()  {
  addVertex(_theZPQQVertex);
  StandardModel::doinit();
}

ZprimeModel::ZprimeModel():_gZPUU_L(1.0), _gZPUU_R(1.0), _gZPDD_L(1.0), _gZPDD_R(1.0), _gZPCC_L(1.0), _gZPCC_R(1.0), _gZPSS_L(1.0),  _gZPSS_R(1.0), ,_gZPBB_L(1.0), _gZPBB_R(1.0),_gZPTT_L(1.0), _gZPTT_R(1.0),_ZPoverall(1.0) {}

IBPtr ZprimeModel::clone() const {
  return new_ptr(*this);
}
IBPtr ZprimeModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ZprimeModel::persistentOutput(PersistentOStream & os) const {
  os  << _theZPQQVertex    
     << _gZPTT_L
     << _gZPTT_R
     << _gZPBB_L
     << _gZPBB_R
     << _gZPCC_L
     << _gZPCC_R
     << _gZPSS_L
     << _gZPSS_R
     << _gZPUU_L
     << _gZPUU_R
     << _gZPDD_L
     << _gZPDD_R
     << _ZPoverall;
}

void ZprimeModel::persistentInput(PersistentIStream & is, int) {
  is  >> _theZPQQVertex    
     >> _gZPTT_L
     >> _gZPTT_R
     >> _gZPBB_L
     >> _gZPBB_R
     >> _gZPCC_L
     >> _gZPCC_R
     >> _gZPSS_L
     >> _gZPSS_R
     >> _gZPUU_L
     >> _gZPUU_R
     >> _gZPDD_L
     >> _gZPDD_R
     >> _ZPoverall;
}

ClassDescription<ZprimeModel> ZprimeModel::initZprimeModel;
// Definition of the static class description member.

void ZprimeModel::Init() {

 static Reference<ZprimeModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexZPQQ
  ("Vertex/ZPQQ",
   "Reference to the Z prime Quark-Antiquark vertex",
   &ZprimeModel::_theZPQQVertex, false, false, true, false, false);



  static Parameter<ZprimeModel, double> interfaceZPTTCoupling
    ("ZPTTLCoupling",
     "The left-handed Z prime coupling to top anti-top",
     &ZprimeModel::_gZPTT_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPTTRCoupling
    ("ZPTTRCoupling",
     "The right-handed Z prime coupling to top anti-top",
     &ZprimeModel::_gZPTT_R, 1.0, 0., 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPUULCoupling
    ("ZPUULCoupling",
     "The left-handed Z prime coupling to up upbar",
     &ZprimeModel::_gZPUU_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPUURCoupling
    ("ZPUURCoupling",
     "The right-handed Z prime coupling to up upbar",
     &ZprimeModel::_gZPUU_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPCCLCoupling
    ("ZPCCLCoupling",
     "The left-handed Z prime coupling to c cbar",
     &ZprimeModel::_gZPCC_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPCCRCoupling
    ("ZPCCRCoupling",
     "The right-handed Z prime coupling to c cbar",
     &ZprimeModel::_gZPCC_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  
  static Parameter<ZprimeModel, double> interfaceZPBBCoupling
    ("ZPBBLCoupling",
     "The left-handed Z prime coupling to b bbar",
     &ZprimeModel::_gZPBB_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPBBRCoupling
    ("ZPBBRCoupling",
     "The right-handed Z prime coupling to b bbar",
     &ZprimeModel::_gZPBB_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);


  static Parameter<ZprimeModel, double> interfaceZPDDLCoupling
    ("ZPDDLCoupling",
     "The left-handed Z prime coupling to d dbar",
     &ZprimeModel::_gZPDD_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPDDRCoupling
    ("ZPDDRCoupling",
     "The right-handed Z prime coupling to d dbar",
     &ZprimeModel::_gZPDD_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

 static Parameter<ZprimeModel, double> interfaceZPSSLCoupling
    ("ZPSSLCoupling",
     "The left-handed Z prime coupling to s sbar",
     &ZprimeModel::_gZPSS_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ZprimeModel, double> interfaceZPSSRCoupling
    ("ZPSSRCoupling",
     "The right-handed Z prime coupling to s sbar",
     &ZprimeModel::_gZPSS_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
 static Parameter<ZprimeModel, double> interfaceZPOverallCoupling
    ("ZPOverallCoupling",
     "An overall coupling of the Z prime to quark-antiquark",
     &ZprimeModel::_ZPoverall, 1.0, -10000.0, 10000.0,
     false, false, Interface::limited);
  
  static ClassDocumentation<ZprimeModel> documentation
    ("There is no documentation for the ZprimeModel class");

}

