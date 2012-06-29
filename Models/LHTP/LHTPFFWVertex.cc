// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFWVertex class.
//

#include "LHTPFFWVertex.h"
#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/CKMBase.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

IBPtr LHTPFFWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFWVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _ckm << _salpha << _calpha << _sbeta << _cbeta;
}

void LHTPFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _ckm >> _salpha >> _calpha >> _sbeta >> _cbeta;
}

ClassDescription<LHTPFFWVertex> 
LHTPFFWVertex::initLHTPFFWVertex;
// Definition of the static class description member.

void LHTPFFWVertex::Init() {

  static ClassDocumentation<LHTPFFWVertex> documentation
    ("The LHTPFFWVertex class implements the couplings of the W"
     " and W_H bosons of the Little Higgss model with T-parity to the fermions.");

}

void LHTPFFWVertex::doinit() {
  ThePEG::Helicity::FFVVertex::doinit();
  Ptr<CKMBase>::transient_pointer CKM = generator()->standardModel()->CKM();
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(CKM);
  if(hwCKM) {
    vector< vector<Complex > > CKM;
    CKM = hwCKM->getUnsquaredMatrix(generator()->standardModel()->families());
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int iy=0;iy<3;++iy) {
	_ckm[ix][iy]=CKM[ix][iy];
      }
    }
  }
  else {
    throw InitException() << "Must have access to the Herwig::StandardCKM object"
			  << "for the CKM matrix in LHTPFFWVertex::doinit()"
			  << Exception::runerror;
  }
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFWVertex::doinit()"
				   << Exception::runerror;
  _salpha = model->sinAlpha();
  _calpha = model->cosAlpha();
  _sbeta  = model->sinBeta();
  _cbeta  = model->cosBeta();
  //  cerr << *this;
}

LHTPFFWVertex::LHTPFFWVertex()
  : _ckm(3,vector<Complex>(3,0.0)), _couplast(0.),_q2last(0.*sqr(MeV))  {
  orderInGem(1);
  orderInGs(0);
  // particles for the vertex
  vector<long> first,second,third;
  // particles for outgoing W-
  // quarks
  for(int ix = 1; ix < 6; ix += 2) {
    for(int iy = 2; iy < 7; iy += 2) {
      addToList(-ix,      iy,      -24);
    }
  }
  // additional T quark
  addToList(-5,   8,  -24);
  // leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix,    ix + 1,    -24);
  }
  // T-odd quarks
  for(long ix = 4000002; ix < 4000007; ix += 2) {
    addToList(-ix + 1,    ix,    -24);
  }
  // T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix,    ix + 1,    -24);
  }
  // T-odd leptons
  // particles for outgoing W+
  // quarks
  for(int ix = 2; ix < 7; ix += 2) {
    for(int iy = 1; iy < 6; iy += 2) {
      addToList(-ix,       iy,       24);
    }
  }
  // additional T quark
  addToList(-8,   5,  24);
  // leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix - 1,    ix,    24);
  }
  // T-odd quarks
  for(long ix = 4000002; ix < 4000009; ix += 2) {
    addToList(-ix,    ix-1,    24);
  }
  // T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix-1,    ix,    24);
  }
  // particles for w_H-
  // quark and T-odd quark
  for(int ix = 1; ix < 6; ix += 2) {
    addToList(-ix-4000000,    ix+1,    -34);
    addToList(-ix,    ix+1+4000000,    -34);
  }
  // lepton and T-odd lepton
  for(int ix = 11;ix < 17; ix += 2) {
    addToList(-ix-4000000,    ix+1,    -34);
    addToList(-ix,    ix+1+4000000,    -34);
  }
  // particles for w_h+
  // quark and T-odd quark
  for(int ix = 1;ix < 6;ix += 2) {
    addToList(ix + 4000000,    -ix - 1,    34);
    addToList(ix,    -ix - 4000001,     34);
  }
  // leptons and T-odd lepton
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix - 4000001,    ix,    34);
    addToList(-ix - 1,    ix + 4000000,    34);
  }
}

// coupling for FFW vertex
void LHTPFFWVertex::setCoupling(Energy2 q2, tcPDPtr a,
				tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -sqrt(0.5)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // SM W boson
  if(abs(c->id())==ParticleID::Wplus) {
    long ia(abs(a->id())),ib(abs(b->id()));
    // quarks
    if(ia >= 1 && ia <= 6 && ib >= 1 && ib <= 6 ) {
      int iu,id;
      // up type first
      if(ia % 2 == 0) {
	iu = ia/2;
	id = (ib+1)/2;
      }
      // down type first
      else {
	iu = ib/2;
	id = (ia+1)/2;
      }
      if( iu<1 || iu>3 || id<1 || id>3)
	throw HelicityConsistencyError() << "LHTPFFWVertex::setCoupling "
					 << "Unknown particle in W vertex" 
					 << Exception::runerror;
      if(ia != 6 && ib != 6) left(_ckm[iu-1][id-1]);
      else left(_ckm[iu-1][id-1]*_cbeta);
      right(0.);
    }
    else if(ia == 8|| ib == 8) {
      left(-_sbeta);
      right(0.);
    }
    // leptons
    else if( ia >= 11 && ia <= 16) {
      left(1.);
      right(0.);
    }
    else {
      left(1.);
      right(1.);
    }
  }
  else {
    left(1.);
    right(0.);
  }
}
