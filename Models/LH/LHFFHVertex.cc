// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFHVertex class.
//

#include "LHFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coup,1./GeV) << _model;
}

void LHFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,1./GeV) >> _model;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFHVertex,FFSVertex>
describeHerwigLHFFHVertex("Herwig::LHFFHVertex", "HwLHModel.so");

void LHFFHVertex::Init() {

  static ClassDocumentation<LHFFHVertex> documentation
    ("The LHFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model");

}

LHFFHVertex::LHFFHVertex() 
  : _q2last(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  _masslast[0] = 0.*GeV; 
  _masslast[1] = 0.*GeV;
  _idlast[0] = 0;
  _idlast[1] = 0;
  colourStructure(ColourStructure::DELTA);
}

void LHFFHVertex::doinit() {
  // SM like higgs
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  25);
  }
  addToList(  -6,   8,  25);
  addToList(  -8,   6,  25);
  addToList(  -8,   8,  25);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  25);
  }
  // phi0
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  35);
  }
  addToList(  -6,   8,  35);
  addToList(  -8,   6,  35);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  35);
  }
  // phiP
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  36);
  }
  addToList(  -6,   8,  36);
  addToList(  -8,   6,  36);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  36);
  }
  // phi +/-
  for(int ix=1;ix<6;ix+=2) {
    addToList( -ix-1,   ix,  37);
    addToList( -ix  , ix+1, -37);
  }
  addToList( -8 ,   5,  37);
  addToList( -5 ,   8, -37);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix-1,   ix,  37);
    addToList( -ix  , ix+1, -37);
  }
  _model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!_model)   throw InitException() << "Must be using the LHModel "
				      << " in LHFFPVertex::doinit()"
				      << Exception::runerror;
  _coup.resize(11);
  Energy v     = _model->vev();
  double s0    = _model->sinTheta0();
  double sP    = _model->sinThetaP();
  double sPlus = _model->sinThetaPlus();
  double s02   = sqr(s0);
  double vf    = _model->vev()/_model->f();
  double xL    = sqr(_model->lambda1())/(sqr(_model->lambda1())+sqr(_model->lambda2()));
  double xR    = sqr(_model->lambda1())/sqrt(sqr(_model->lambda1())+sqr(_model->lambda2()));
  Energy mT    = getParticleData(8)->mass();
  // lightest higgs couplings
  // coupling of light SM fermions
  _coup[0] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf))/v;
  // couplings to top quark
  _coup[1] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf)+sqr(vf)*xL*(1.+xL))/v;
  // couplings to the T quark
  _coup[2] =-xR*(1.+xL)*vf/mT;
  // couplings to tT
  _coup[3] = xR/mT;
  _coup[4] = vf/v*(1.+xL);
  // phi 0
  // light particles
  _coup[5] = sqrt(0.5)/v*(vf-sqrt(2.)*s0);
  // mixed
  _coup[6] = sqrt(0.5)/v*(vf-sqrt(2.)*s0)*_model->lambda1()/_model->lambda2();
  // phi P
  _coup[7] = Complex(0.,1.)*sqrt(0.5)/v*(vf-sqrt(2.)*sP);
  _coup[8] = Complex(0.,1.)*sqrt(0.5)/v*(vf-sqrt(2.)*sP)*_model->lambda1()/_model->lambda2();
  // phi +/-
  _coup[9] = -sqrt(0.5)/v*(vf-2.*sPlus);
  _coup[9] = -sqrt(0.5)/v*(vf-2.*sPlus)*_model->lambda1()/_model->lambda2();
  FFSVertex::doinit();
}

void LHFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  // left and right couplings set to one
  left (1.);
  right(1.);
  // first the overall normalisation
  if(q2!=_q2last||_idlast[0]!=iferm||_idlast[1]!=ianti) {
    _q2last = q2;
    _idlast[0] = iferm;
    if(_idlast[0]==8) _idlast[0]=6;
    assert((_idlast[0]>=1  && _idlast[0]<=6 ) || 
	   (_idlast[0]>=11 && _idlast[0]<=16));
    if(iferm==_idlast[0])
      _masslast[0] = _model->mass(q2,a);
    else
      _masslast[0] = _model->mass(q2,getParticleData(ParticleID::t));
    _idlast[1] = ianti;
    if(_idlast[1]==8) _idlast[1]=6;
    assert((_idlast[1]>=1  && _idlast[1]<=6 ) || 
	   (_idlast[1]>=11 && _idlast[1]<=16));
    if(_idlast[0]!=_idlast[1]) {
      if(ianti==_idlast[1])
	_masslast[1] = _model->mass(q2,a);
      else
	_masslast[1] = _model->mass(q2,getParticleData(ParticleID::t));
    }
    else {
      _masslast[1] = _masslast[0];
    }
  }
  // SM like higgs
  if(c->id()==ParticleID::h0) {
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=5 ) || 
	 (iferm>=11 && iferm<=16)) {
	norm(-Complex(_coup[0]*_masslast[0]));
      }
      else if(iferm==6) {
	norm(-Complex(_coup[1]*_masslast[0]));
      }
      else if(iferm==8) {
	norm(-Complex(_coup[2]*a->mass()));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (ianti == 6 && iferm == 8 ));
      Complex cleft,cright;
      if(iferm==6) {
	cleft  = Complex(-_coup[3]*b->mass());
	cright = Complex(-_coup[4]*_masslast[0]);
      }
      else {
	cleft  = Complex(-_coup[3]*a->mass());
	cright = Complex(-_coup[4]*_masslast[0]);
      }
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::H0) {
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=6 ) || 
	 (iferm>=11 && iferm<=16)) {
	norm(-Complex(_coup[5]*_masslast[0]));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (iferm == 8 && ianti == 6 ) );
      Complex cleft  = Complex(_coup[6]*_masslast[0]);
      Complex cright = 0.;
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::A0) {
    left(-1.);
    right(1.);
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=6 ) || 
	 (iferm>=11 && iferm<=16)) {
	if(iferm%2==0)
	  norm(-Complex( _coup[7]*_masslast[0]));
	else
	  norm(-Complex(-_coup[7]*_masslast[0]));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (iferm == 8 && ianti == 6 ));
      Complex cleft  = Complex(_coup[8]*_masslast[0]);
      Complex cright = 0.;
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::Hplus) {
    norm(1.);
    Complex cleft(0.),cright(0.);
    if(iferm%2==0) {
      if(iferm==ParticleID::t) {
	cleft  = Complex(_masslast[0]*_coup[ 9]);
      }
      else {
	cleft  = Complex(_masslast[0]*_coup[10]);
	cright = Complex(_masslast[1]*_coup[10]);
      }
    }
    else {
      if(ianti==ParticleID::t) {
	cleft  = Complex(_masslast[1]*_coup[ 9]);
      }
      else {
	cleft  = Complex(_masslast[1]*_coup[10]);
	cright = Complex(_masslast[0]*_coup[10]);
      }
    }
    left ( cleft);
    right(cright);
  }
  else if(c->id()==ParticleID::Hminus) {
    norm(1.);
    Complex cleft(0.),cright(0.);
    if(iferm%2==0) {
      if(iferm==ParticleID::t) {
	cright = Complex(_masslast[0]*_coup[ 9]);
      }
      else {
	cright = Complex(_masslast[0]*_coup[10]);
	cleft  = Complex(_masslast[1]*_coup[10]);
      }
    }
    else {
      if(ianti==ParticleID::t) {
	cright = Complex(_masslast[1]*_coup[ 9]);
      }
      else {
	cright = Complex(_masslast[1]*_coup[10]);
	cleft  = Complex(_masslast[0]*_coup[10]);
      }
    }
    left ( cleft);
    right(cright);
  }
}
