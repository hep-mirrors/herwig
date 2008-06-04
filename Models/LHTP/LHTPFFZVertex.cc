// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFZVertex class.
//

#include "LHTPFFZVertex.h"
#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHTPFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr << _tl << _tr << _coupd << _coupu << _coupe << _coupnu ;
}

void LHTPFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _tl >> _tr >> _coupd >> _coupu >> _coupe >> _coupnu;
}

ClassDescription<LHTPFFZVertex> LHTPFFZVertex::initLHTPFFZVertex;
// Definition of the static class description member.

void LHTPFFZVertex::Init() {

  static ClassDocumentation<LHTPFFZVertex> documentation
    ("There is no documentation for the LHTPFFZVertex class");

}

LHTPFFZVertex::LHTPFFZVertex() : _gl(37,0.0), _gr(37,0.0),
				 _tl(5 ,0.0), _tr(5 ,0.0),
				 _coupd(0.), _coupu(0.),
				 _coupe(0.), _coupnu(0.),
				 _couplast(0.0), _q2last(0.*GeV2) {
  // PDG codes for the particles
  vector<long> first,second,third;
  // Z
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
  }
  // T+T+
  first .push_back(-8);
  second.push_back(+8);
  third.push_back(23);
  //T+t
  first .push_back(-6);
  second.push_back(+8);
  third.push_back(23);
  first .push_back(-8);
  second.push_back(+6);
  third.push_back(23);
  // the leptons
  for(int ix=11;ix<17;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
  }
  // the T-odd quarks
  for(long ix=1;ix<7;++ix) {
    first.push_back(-ix-4000000);
    second.push_back(ix+4000000);
    third.push_back(23);
  }
  first .push_back(-4000008);
  second.push_back(+4000008);
  third.push_back(23);
  // the T-odd leptons
  for(long ix=11;ix<17;++ix) {
    first.push_back(-ix-4000000);
    second.push_back(ix+4000000);
    third.push_back(23);
  }
  // Z_H
  // the quarks
  for(long ix=1;ix<7;++ix) {
    first.push_back(-ix-4000000);
    second.push_back(ix);
    third.push_back(33);
    first.push_back(-ix);
    second.push_back(ix+4000000);
    third.push_back(33);
  }  
  first .push_back(      -8);
  second.push_back(+4000008);
  third .push_back(33);
  first .push_back(       8);
  second.push_back(-4000008);
  third .push_back(33);
  first .push_back(      -6);
  second.push_back(+4000008);
  third .push_back(33);
  first .push_back(       6);
  second.push_back(-4000008);
  third .push_back(33);
  // the leptons
  for(long ix=11;ix<17;++ix) {
    first.push_back(-ix-4000000);
    second.push_back(ix);
    third.push_back(33);
    first.push_back(-ix);
    second.push_back(ix+4000000);
    third.push_back(33);
  }
  setList(first,second,third);
}

void LHTPFFZVertex::doinit() throw(InitException) {
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
  orderInGem(1);
  orderInGs(0);
  FFVVertex::doinit();
}

void LHTPFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm = abs(a->id());
  int ianti = abs(b->id());
  int ibos  = c->id();
  if(ibos==ParticleID::Z0) {
    if(iferm==8||ianti==8) {
      if(iferm==8&&ianti==8) {
	setLeft (_tl[1]);
	setRight(_tr[1]);
      }
      else {
	setLeft (_tl[2]);
	setRight(_tr[2]);
      }
    }
    else if(iferm==4000008) {
      setLeft (_tl[0]);
      setRight(_tr[0]);
    }
    else if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      setLeft(_gl[iferm]);
      setRight(_gr[iferm]);
    }
    else {
      iferm = iferm%4000000+20;
      setLeft (_gl[iferm]);
      setRight(_gr[iferm]);
    }
  }
  else {
    if(iferm==4000008||ianti==4000008) {
      if(iferm==8||ianti==8) {
	setLeft (_tl[4]);
	setRight(_tr[4]);
      }
      else {
	setLeft (_tl[3]);
	setRight(_tr[3]);
      }
    }
    else {
      iferm = iferm%4000000;
      setRight(0.);
      if(iferm<=6) {
	if(iferm%2==0) setLeft(_coupu );
	else           setLeft(_coupd );
      }
      else {
	if(iferm%2==0) setLeft(_coupnu);
	else           setLeft(_coupe );
      }
    }
  }
}
