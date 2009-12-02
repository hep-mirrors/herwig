// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaZPrimeZVertex class.
//

#include "GammaZPrimeZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RadiativeZPrimeModel.h"

using namespace RadiativeZPrime;

void GammaZPrimeZVertex::persistentOutput(PersistentOStream & os) const {
  os << _coup;
}

void GammaZPrimeZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _coup;
}

ClassDescription<GammaZPrimeZVertex> GammaZPrimeZVertex::initGammaZPrimeZVertex;
// Definition of the static class description member.

void GammaZPrimeZVertex::Init() {

  static ClassDocumentation<GammaZPrimeZVertex> documentation
    ("The GammaZPrimeZVertex class implements the anomalous Z Z' gamma"
     " vertex in the radiativeZPrimeModel");

}

void GammaZPrimeZVertex::setCoupling(Energy2, tcPDPtr ,tcPDPtr,tcPDPtr) {
  norm(_coup);
}

GammaZPrimeZVertex::GammaZPrimeZVertex() {
  // PDG codes for the particles
  addToList(22,32,23);
}

void GammaZPrimeZVertex::doinit() {
  tcSMPtr sm = generator()->standardModel();
  tcRadiativeZPrimeModelPtr model = 
    dynamic_ptr_cast<tcRadiativeZPrimeModelPtr>(generator()->standardModel());
  // calculate the coupling
  double F=0.5*(2./3.+1./3.);
  Energy2 scale = sqr(getParticleData(32)->mass());
  double sw = sqrt(model->sin2ThetaW());
  double cw = sqrt(1.-model->sin2ThetaW());
  _coup = 1./16/sqr(Constants::pi)*4.*Constants::pi*model->alphaEM(scale)/sw/cw*
    model->gZprime()*F;
  orderInGem(3);
  orderInGs(0);
  AbstractVVVVertex::doinit();
}
