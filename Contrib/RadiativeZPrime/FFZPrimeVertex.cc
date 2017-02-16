// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFZPrimeVertex class.
//

#include "FFZPrimeVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RadiativeZPrimeModel.h"

using namespace RadiativeZPrime;

void FFZPrimeVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr;
}

void FFZPrimeVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr;
}

ClassDescription<FFZPrimeVertex> FFZPrimeVertex::initFFZPrimeVertex;
// Definition of the static class description member.

void FFZPrimeVertex::Init() {

  static ClassDocumentation<FFZPrimeVertex> documentation
    ("There is no documentation for the FFZPrimeVertex class");

}

void FFZPrimeVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  norm(1.);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    left(_gl[iferm]);
    right(_gr[iferm]);
  }
  else
    throw HelicityConsistencyError() << "FFZPrimeVertex::setCoupling "
				     << "Unknown particle in Z vertex" 
				     << Exception::runerror;
}

FFZPrimeVertex::FFZPrimeVertex() : _gl(17,0.0), _gr(17,0.0) {
  // PDG codes for the particles
  // the quarks
  for(long ix=1;ix<7;++ix) {
    addToList(-ix,ix,32);
  }
  // the leptons
  for(long ix=11;ix<17;++ix) {
    addToList(-ix,ix,32);
  }
}

void FFZPrimeVertex::doinit() {
  tcSMPtr sm = generator()->standardModel();
  tcRadiativeZPrimeModelPtr model = 
    dynamic_ptr_cast<tcRadiativeZPrimeModelPtr>(generator()->standardModel());
  double fact = 0.25*model->gZprime();
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = fact*(model->zPrimevd()  + model->zPrimead() );
    _gl[2*ix ]   = fact*(model->zPrimevu()  + model->zPrimeau() );
    _gl[2*ix+9 ] = fact*(model->zPrimeve()  + model->zPrimeae() );
    _gl[2*ix+10] = fact*(model->zPrimevnu() + model->zPrimeanu());
    _gr[2*ix-1]  = fact*(model->zPrimevd()  - model->zPrimead() );
    _gr[2*ix ]   = fact*(model->zPrimevu()  - model->zPrimeau() );
    _gr[2*ix+9 ] = fact*(model->zPrimeve()  - model->zPrimeae() );
    _gr[2*ix+10] = fact*(model->zPrimevnu() - model->zPrimeanu());
  }
  orderInGem(1);
  orderInGs(0);
  FFVVertex::doinit();
}
