// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsFFZVertex class.
//

#include "LittleHiggsFFZVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LittleHiggsFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _model << _gl << _gr << _glH << _grH;
}

void LittleHiggsFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _gl >> _gr >> _glH >> _grH;
}

ClassDescription<LittleHiggsFFZVertex> LittleHiggsFFZVertex::initLittleHiggsFFZVertex;
// Definition of the static class description member.

void LittleHiggsFFZVertex::Init() {

  static ClassDocumentation<LittleHiggsFFZVertex> documentation
    ("There is no documentation for the LittleHiggsFFZVertex class");

}

LittleHiggsFFZVertex::LittleHiggsFFZVertex() : _couplast(0.0), _q2last(0.*GeV2) {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(33);
  }
  first.push_back( -8);
  second.push_back( 8);
  third.push_back(23);
  first.push_back( -8);
  second.push_back( 8);
  third.push_back(33);
  first.push_back( -6);
  second.push_back( 8);
  third.push_back(23);
  first.push_back( -8);
  second.push_back( 6);
  third.push_back(33);
  // the leptons
  for(unsigned int ix=11;ix<17;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(33);
  }
  setList(first,second,third);
}

void LittleHiggsFFZVertex::doinit() throw(InitException) {
  FFVVertex::doinit();
  _model = dynamic_ptr_cast<cLittleHiggsModelPtr>(generator()->standardModel());
  if(!_model) 
    throw InitException() << "Must be using the LittleHiggsModel "
			  << " in LittleHiggsFFZVertex::doinit()"
			  << Exception::runerror;
  double sw2(_model->sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double pre =-0.5/sqrt(sw2*(1.-sw2));
  double s (_model->sinTheta()     ),c (_model->cosTheta()     );
  double sp(_model->sinThetaPrime()),cp(_model->cosThetaPrime());
  double sp2(sqr(sp)),cp2(sqr(cp));
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  double xB(-2.5/sw*sp*cp*(sqr(cp)-sqr(sp)));
  double yu  = -0.4, ye  =  0.6;
  double vf(_model->vev()/_model->f());
  double xL(sqr(_model->lambda1())/(sqr(_model->lambda1())+sqr(_model->lambda2())));
  double vu  = pre*( 0.5-4./3.*sw2-sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+7./15.  -cp2/6.))); 
  double vd  = pre*(-0.5+2./3.*sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+11./15. +cp2/6.))); 
  double ve  = pre*(-0.5+2.*   sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*ye-9./5. +1.5*cp2))); 
  double vv  = pre*(+0.5          -sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(   ye-4./5. +0.5*cp2))); 
  double au  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(+0.2-0.5*cp2))); 
  double ad  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s+sw*xB/sp/cp*(-0.2+0.5*cp2))); 
  double ae  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s+sw*xB/sp/cp*(-0.2+0.5*cp2))); 
  double av  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(+0.2-0.5*cp2))); 
  
  double vtl  = pre*( 0.5-4./3.*sw2-sqr(vf)*(-0.5*sqr(xL)+0.5*cw*xW*c/s
					     +sw*xB/sp/cp*(2.*yu+9./5.-1.5*cp2
							   +(7./15.-2.*cp2/3.)*xL)));
  double atl  = pre*(-0.5-sqr(vf)*(+0.5*sqr(xL)-0.5*cw*xW*c/s
				   +sw*xB/sp/cp*(+0.2-0.5*cp2-0.2*xL)));
  double vth = 2./3.*sw/cw;
  double ath = 0.;
  double vtm = 0.25*xL*vf/cw/sw;
  double atm = -vtm;
  _gl.resize(17);
  _gr.resize(17);
  for(unsigned ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = vd - ad;
    _gl[2*ix ]   = vu - au;
    _gl[2*ix+9 ] = ve - ae;
    _gl[2*ix+10] = vv - av;
    _gr[2*ix-1]  = vd + ad;
    _gr[2*ix ]   = vu + au;
    _gr[2*ix+9 ] = ve + ae;
    _gr[2*ix+10] = vv + av;
  }
  _gl[6] = vtl + atl;
  _gr[6] = vtl - atl;
  _gl[7] = vtm + atm;
  _gr[7] = vtm - atm;
  _gl[8] = vth + ath;
  _gr[8] = vth - ath;
  // heavy Z
  vu  =  0.25*c/s/sw;
  vd  = -0.25*c/s/sw;
  ve  = -0.25*c/s/sw;
  vv  =  0.25*c/s/sw;
  au  = -0.25*c/s/sw;
  ad  =  0.25*c/s/sw;
  ae  =  0.25*c/s/sw;
  av  = -0.25*c/s/sw;
  vtl =  0.25*c/s/sw;
  atl = -0.25*c/s/sw;
  vth =  0.;
  ath =  0.;
  vtm =  -0.25*xL*vf*c/s/sw;
  atm = -vtm;
  _glH.resize(17);
  _grH.resize(17);
  for(unsigned ix=1;ix<4;++ix) {
    _glH[2*ix-1]  = vd - ad;
    _glH[2*ix ]   = vu - au;
    _glH[2*ix+9 ] = ve - ae;
    _glH[2*ix+10] = vv - av;
    _grH[2*ix-1]  = vd + ad;
    _grH[2*ix ]   = vu + au;
    _grH[2*ix+9 ] = ve + ae;
    _grH[2*ix+10] = vv + av;
  }
  _glH[6] = vtl + atl;
  _grH[6] = vtl - atl;
  _glH[7] = vtm + atm;
  _grH[7] = vtm - atm;
  _glH[8] = vth + ath;
  _grH[8] = vth - ath;
  // set order in the couplings
  orderInGem(1);
  orderInGs(0);
}

void LittleHiggsFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _model->alphaEM(q2);
    _couplast = -sqrt(4.0*Constants::pi*alpha);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm = abs(a->id());
  int ianti = abs(b->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)|| iferm == 8) {
    // Z0
    if(c->id()==ParticleID::Z0) {
      if(ianti==iferm) {
	setLeft(_gl[iferm]);
	setRight(_gr[iferm]);
      }
      else {
	setLeft (_gl[7]);
	setRight(_gr[7]);
      }
    }
    else {
      if(ianti==iferm) {
	setLeft (_glH[iferm]);
	setRight(_grH[iferm]);
      }
      else {
	setLeft (_glH[7]);
	setRight(_grH[7]);
      }
    }
  }
  else
    throw HelicityConsistencyError() << "LittleHiggsFFZVertex::setCoupling "
				     << "Unknown particle in Z vertex"
				     << a->PDGName() << " " << b->PDGName()
				     << " " << c->PDGName()
				     << Exception::runerror;
}
