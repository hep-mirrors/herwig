// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFZVertex class.
//

#include "LHTPFFZVertex.h"
#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPFFZVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFZVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr << _tl << _tr << _coupd << _coupu << _coupe << _coupnu ;
}

void LHTPFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _tl >> _tr >> _coupd >> _coupu >> _coupe >> _coupnu;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFZVertex,FFVVertex>
describeHerwigLHTPFFZVertex("Herwig::LHTPFFZVertex", "HwLHTPModel.so");

void LHTPFFZVertex::Init() {

  static ClassDocumentation<LHTPFFZVertex> documentation
    ("The LHTPFFZVertex class implements the couplings of "
     "the fermions to the Z boson and its heavy partner in the"
     " Little Higgs model with T-parity.");

}

LHTPFFZVertex::LHTPFFZVertex() : _gl(37,0.0), _gr(37,0.0),
				 _tl(5 ,0.0), _tr(5 ,0.0),
				 _coupd(0.), _coupu(0.),
				 _coupe(0.), _coupnu(0.),
				 _couplast(0.0), _q2last(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  // Z
  // the quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix,    23);
  }
  // T+T+
  addToList(-8,  +8,  23);
  //T+t
  addToList(-6,  +8,  23);
  addToList(-8,  +6,  23);
  // the leptons
  for(int ix = 11; ix < 17; ++ix) {
    addToList(-ix,    ix,    23);
  }
  // the T-odd quarks
  for(long ix = 4000001; ix < 4000007; ++ix) {
    addToList(-ix,    ix,    23);
  }
  addToList(-4000008,  +4000008,  23);
  // the T-odd leptons
  for(long ix = 4000011;ix<17;++ix) {
    addToList(-ix-4000000,    ix+4000000,    23);
  }
  // Z_H
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }  
  addToList(      -8,  +4000008,  33);
  addToList(       8,  -4000008,  33);
  addToList(      -6,  +4000008,  33);
  addToList(       6,  -4000008,  33);
  // the leptons
  for(unsigned int ix=11;ix<17;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }
}

void LHTPFFZVertex::doinit() {
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFPVertex::doinit()"
				   << Exception::runerror;
  double sw = model->sin2ThetaW();
  double cw = sqrt(1.-sw);
  sw = sqrt(sw);
  double fact = 0.25/sw/cw;
  for(int ix=1;ix<4;++ix) {
    // SM fermions
    _gl[2*ix-1]  = fact*(model->vd()  + model->ad() );
    _gl[2*ix ]   = fact*(model->vu()  + model->au() );
    _gl[2*ix+9 ] = fact*(model->ve()  + model->ae() );
    _gl[2*ix+10] = fact*(model->vnu() + model->anu());
    _gr[2*ix-1]  = fact*(model->vd()  - model->ad() );
    _gr[2*ix ]   = fact*(model->vu()  - model->au() );
    _gr[2*ix+9 ] = fact*(model->ve()  - model->ae() );
    _gr[2*ix+10] = fact*(model->vnu() - model->anu());
    // T-odd fermions
    _gl[2*ix-1 +20] = fact*(model->vd()  + model->ad() );
    _gl[2*ix   +20] = fact*(model->vu()  + model->au() );
    _gl[2*ix+9 +20] = fact*(model->ve()  + model->ae() );
    _gl[2*ix+10+20] = fact*(model->vnu() + model->anu());
    _gr[2*ix-1 +20] = _gl[2*ix-1 +20];
    _gr[2*ix   +20] = _gl[2*ix   +20];
    _gr[2*ix+9 +20] = _gl[2*ix+9 +20];
    _gr[2*ix+10+20] = _gl[2*ix+10+20];
  }
  // couplngis to Z for extended top sector
  _tl[0] = -2./3.*sw/cw;
  _tr[0] = -2./3.*sw/cw;
  _tl[1] = (0.5*sqr(model->sinBeta())-2./3.*sqr(sw))/cw/sw;
  _tr[1] = -2./3.*sw/cw;
  _tl[2] = -0.5/sw/cw*model->sinBeta()*model->cosBeta();
  _tr[2] = 0.;
  // couplings of fermion T-odd fermion Z_H
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  _coupd  =  0.1*(sH/cw+cH/sw);
  _coupu  =  0.1*(sH/cw-cH/sw);
  _coupe  =  0.1*(sH/cw+cH/sw);
  _coupnu =  0.1*(sH/cw-cH/sw);
  // couplings of T-odd top
  _tl[3] = 0.4*sH*model->sinBeta ()/cw;
  _tr[3] = 0.4*sH*model->sinAlpha()/cw;
  _tl[4] = 0.4*sH*model->cosBeta ()/cw;
  _tr[4] = 0.4*sH*model->cosAlpha()/cw;
  FFVVertex::doinit();
}

void LHTPFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  long iferm = abs(a->id());
  long ianti = abs(b->id());
  long ibos  = c->id();
  if(ibos == ParticleID::Z0) {
    if(iferm == 8 || ianti == 8) {
      if(iferm == 8 && ianti == 8) {
	left (_tl[1]);
	right(_tr[1]);
      }
      else {
	left (_tl[2]);
	right(_tr[2]);
      }
    }
    else if(iferm == 4000008) {
      left (_tl[0]);
      right(_tr[0]);
    }
    else if((iferm >= 1 && iferm <= 6)|| (iferm >= 11 && iferm <= 16)) {
      left(_gl[iferm]);
      right(_gr[iferm]);
    }
    else {
      iferm = (iferm % 4000000) + 20;
      left (_gl[iferm]);
      right(_gr[iferm]);
    }
  }
  else {
    if(iferm == 4000008||ianti == 4000008) {
      if(iferm == 8||ianti == 8) {
	left (_tl[4]);
	right(_tr[4]);
      }
      else {
	left (_tl[3]);
	right(_tr[3]);
      }
    }
    else {
      iferm = iferm % 4000000;
      right(0.);
      if(iferm <= 6) {
	if(iferm % 2 == 0) left(_coupu );
	else  left(_coupd );
      }
      else {
	if(iferm%2 == 0) left(_coupnu);
	else left(_coupe );
      }
    }
  }
}
