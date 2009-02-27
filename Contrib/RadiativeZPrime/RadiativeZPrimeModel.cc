// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeZPrimeModel class.
//

#include "RadiativeZPrimeModel.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace RadiativeZPrime;

void RadiativeZPrimeModel::persistentOutput(PersistentOStream & os) const {
  os << _gZprime << _useZcouplings 
     << _vnu << _ve << _vu << _vd 
     << _anu << _ae << _au << _ad
     << _ffZPrimeVertex << _gammaZPrimeZVertex;
}

void RadiativeZPrimeModel::persistentInput(PersistentIStream & is, int) {
  is >> _gZprime >> _useZcouplings 
     >> _vnu >> _ve >> _vu >> _vd 
     >> _anu >> _ae >> _au >> _ad
     >> _ffZPrimeVertex >> _gammaZPrimeZVertex;
}

ClassDescription<RadiativeZPrimeModel> RadiativeZPrimeModel::initRadiativeZPrimeModel;
// Definition of the static class description member.

void RadiativeZPrimeModel::Init() {

  static ClassDocumentation<RadiativeZPrimeModel> documentation
    ("The RadiativeZPrimeModel class");

  static Switch<RadiativeZPrimeModel,bool> interfaceUseZCouplings
    ("UseZCouplings",
     "Set the fermion couplings to the Z prime to be the same as for the Z",
     &RadiativeZPrimeModel::_useZcouplings, true, false, false);
  static SwitchOption interfaceUseZCouplingsYes
    (interfaceUseZCouplings,
     "Yes",
     "Set Z' couplings to Z couplings",
     true);
  static SwitchOption interfaceUseZCouplingsNo
    (interfaceUseZCouplings,
     "No",
     "Z' couplings are inputted",
     false);
  
  static Parameter<RadiativeZPrimeModel,double> interfacegZPrime
    ("gZPrime",
     "The coupling of the Z'",
     &RadiativeZPrimeModel::_gZprime, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeVnu
    ("Zprimev_nu",
     "The vector coupling between a neutrino and a Z'. ",
     &RadiativeZPrimeModel::_vnu, 1.0, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeVe
    ("Zprimev_e",
     "The vector coupling between a charged lepton and a Z'. ",
     &RadiativeZPrimeModel::_ve, -0.072, 0.0, 0.0, false, false, false);
  
  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeVu
    ("Zprimev_u",
     "The vector coupling between an up-type quark and a Z'. ",
     &RadiativeZPrimeModel::_vu, 0.3813, 0.0, 0.0, false, false, false);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeVd
    ("Zprimev_d",
     "The vector coupling between a down-type quark and a Z'. ",
     &RadiativeZPrimeModel::_vd, -0.6907, 0.0, 0.0, false, false, false);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeAnu
    ("Zprimea_nu",
     "The axial coupling between a neutrino and a Z'. ",
     &RadiativeZPrimeModel::_anu, 1.0, 0.0, 0.0, false, false, false);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeAe
    ("Zprimea_e",
     "The axial coupling between a charged lepton and a Z'. ",
     &RadiativeZPrimeModel::_ae, -1.0, 0.0, 0.0, false, false, false);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeAu
    ("Zprimea_u",
     "The axial coupling between an up-type quark and a Z'. ",
     &RadiativeZPrimeModel::_au, 1.0, 0.0, 0.0, false, false, false);

  static Parameter<RadiativeZPrimeModel,double> interfaceZPrimeAd
    ("Zprimea_d",
     "The axial coupling between a down-type quark and a Z'. ",
     &RadiativeZPrimeModel::_ad, -1.0, 0.0, 0.0, false, false, false);

  static Reference<RadiativeZPrimeModel,ThePEG::Helicity::AbstractFFVVertex> 
    interfaceVertexFFZ
    ("Vertex/FFZPrime",
     "Reference to the Standard Model FFZ Vertex",
     &RadiativeZPrimeModel::_ffZPrimeVertex, false, false, true, false);

  static Reference<RadiativeZPrimeModel,ThePEG::Helicity::AbstractVVVVertex>
    interfaceVertexGammaZPrimeZ
    ("Vertex/GammaZPrimeZ",
     "Reference to the gamma-Z'-Z vertex",
     &RadiativeZPrimeModel::_gammaZPrimeZVertex, false, false, true, false, false);

}

void RadiativeZPrimeModel::doinit() {
  if(_useZcouplings) {
    _vnu = vnu();
    _ve  = ve() ;
    _vu  = vu() ;
    _vd  = vd() ;
    _anu = anu();
    _ae  = ae() ;
    _au  = au() ;
    _ad  = ad() ;
  }
  StandardModel::doinit();
  addVertex(_ffZPrimeVertex);
  addVertex(_gammaZPrimeZVertex);
}
