// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWWWVertex class.
//

#include "LHTPWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "LHTPModel.h"

using namespace Herwig;

LHTPWWWVertex::LHTPWWWVertex() : coupLast_(0.), q2Last_(ZERO),
				 couplings_(3 ,0.) {
  orderInGem(1);
  orderInGs(0);
}

void LHTPWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << couplings_;
}

void LHTPWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> couplings_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWWWVertex,VVVVertex>
describeHerwigLHTPWWWVertex("Herwig::LHTPWWWVertex", "HwLHTPModel.so");

void LHTPWWWVertex::Init() {

  static ClassDocumentation<LHTPWWWVertex> documentation
    ("The LHTPWWWVertex class implements the coupling of three "
     "electroweak gauge bosons and their heavy partners in the "
     "Little Higgs model with T-parity.");

}

void LHTPWWWVertex::doinit() {
  //SM interactions
  addToList( 24, -24,  22);
  addToList( 24, -24,  23);
  //LHTP
  //W_H W_H A_L
  addToList( 34, -34,  22);
  //W_H W_H Z_L
  addToList( 34, -34,  23);
  //W_H W_L A_H
  addToList( 34, -24,  32);
  addToList( 24, -34,  32); 
  //W_H W_L Z_H
  addToList( 34, -24,  33);
  addToList( 24, -34,  33);
  VVVVertex::doinit();
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "LHTPWWWVertex::doinit() - Model pointer must be of LHTPModel"
      << "type, cannot continue without this."
      << Exception::abortnow;
  double sw(sqrt(model->sin2ThetaW()));
  double cw(sqrt(1. - model->sin2ThetaW()));
  //W W Z
  couplings_[0] = cw/sw;
  //W_L W_H A_H
  couplings_[1] = model->sinThetaH()/sw;
  //W_L W_H Z_H
  couplings_[2] = 1./sw;

}

void LHTPWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  if(q2 != q2Last_) {
    coupLast_ = electroMagneticCoupling(q2);
    q2Last_ = q2;
  }
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // get the PDG code for the neutral boson
  long boson(0);
  if(!a->charged())      {
    boson = ida;
    ida = 22;
  }
  else if(abs(ida) !=ParticleID::Wplus) {
    ida = ida > 0 ? 24 : -24;
  }
  if(!b->charged()) {
    boson = idb;
    idb = 22;
  }
  else if(abs(idb) !=ParticleID::Wplus) {
    idb = idb > 0 ? 24 : -24;
  }
  if(!c->charged()) {
    boson = idc;
    idc = 22;
  }
  else if(abs(idc) !=ParticleID::Wplus) {
    idc = idc > 0 ? 24 : -24;
  }
  assert( boson ==22 || boson==23 || boson==32 || boson==33);
  // get the prefactor
  double pre(0.);
  switch (boson) {
  case 22:
    pre = 1.;
    break;
  case 23:
    pre = couplings_[0];
    break;
  case 32:
    pre = couplings_[1];
    break;
  case 33:
    pre = couplings_[2];
    break;
  default:
    assert(false);
  };
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || 
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )      norm( coupLast_*pre);
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) ) norm(-coupLast_*pre);
  else 
    throw Helicity::HelicityConsistencyError() 
      << "LHTPWWWVertex::setCoupling - Incorrect particles in LHTPWWWVertex. " 
      << a->id() << " " << b->id() << " " << c->id() << '\n'
      << Exception::runerror;
}
