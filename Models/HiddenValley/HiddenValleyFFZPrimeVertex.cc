// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenValleyFFZPrimeVertex class.
//

#include "HiddenValleyFFZPrimeVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "EnumParticles.h"
#include "HiddenValleyModel.h"

using namespace Herwig;

void HiddenValleyFFZPrimeVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr << _gql << _gqr << _gPrime;
}

void HiddenValleyFFZPrimeVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _gql >> _gqr >> _gPrime;
}

ClassDescription<HiddenValleyFFZPrimeVertex> HiddenValleyFFZPrimeVertex::initHiddenValleyFFZPrimeVertex;
// Definition of the static class description member.

void HiddenValleyFFZPrimeVertex::Init() {

  static ClassDocumentation<HiddenValleyFFZPrimeVertex> documentation
    ("There is no documentation for the HiddenValleyFFZPrimeVertex class");

}

void HiddenValleyFFZPrimeVertex::setCoupling(Energy2, tcPDPtr a,tcPDPtr,tcPDPtr) {
  norm(_gPrime);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    left (_gl[iferm]);
    right(_gr[iferm]);
    return;
  }
  iferm -= HiddenID::darkGluon;
  if(iferm<=int(_gql.size())) {
    left (_gql[iferm]);
    right(_gqr[iferm]);
    return;
  }
  assert(false);
}

HiddenValleyFFZPrimeVertex::HiddenValleyFFZPrimeVertex() : _gl(17,0.0), _gr(17,0.0) {
  orderInGem(1);
  orderInGs(0);
  // PDG codes for the particles
  vector<long int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) 
    addToList(-ix,ix,32);
  // the leptons
  for(unsigned int ix=11;ix<17;++ix) 
    addToList(-ix,ix,32);
}

void HiddenValleyFFZPrimeVertex::doinit() {
  tcHiddenValleyPtr model = 
    dynamic_ptr_cast<tcHiddenValleyPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the HiddenValleyModel in "
			  << "HiddenValleyFFZPrimeVertex::doinit()" 
			  << Exception::runerror;
  _gPrime = model->gPrime();
  // SM couplings
  _gl.resize(17,0.);
  _gr.resize(17,0.);
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = model->qL();
    _gl[2*ix ]   = model->qL();
    _gl[2*ix+9 ] = model->lL();
    _gl[2*ix+10] = model->lL();
    _gr[2*ix-1]  = model->dR();
    _gr[2*ix ]   = model->uR();
    _gr[2*ix+9 ] = model->lR();
    _gr[2*ix+10] = 0.  ;
  }
  FFVVertex::doinit();
  _gql.resize(model->NF()+1,0.);
  _gqr.resize(model->NF()+1,0.);
  for(int ix=0;ix<int(model->NF());++ix) {
    int id = HiddenID::darkGluon+1+ix;
    addToList(-id,id,32);
    _gql[ix+1] = model->qCharge()[ix];
    _gqr[ix+1] = model->qCharge()[ix]-2.;
  }
}
