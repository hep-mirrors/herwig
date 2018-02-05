// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFZVertex class.
//

#include "LHFFZVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr << _glH << _grH;
}

void LHFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _glH >> _grH;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFZVertex,FFVVertex>
describeHerwigLHFFZVertex("Herwig::LHFFZVertex", "HwLHModel.so");

void LHFFZVertex::Init() {

  static ClassDocumentation<LHFFZVertex> documentation
    ("The LHFFZVertex class implements the couplings of the Z and Z_H in"
     " the Little Higgs model to the fermions, both of the Standard Model"
     " and the additional heavy top.");

}

LHFFZVertex::LHFFZVertex() : _couplast(0.0), _q2last(0.*GeV2) {
  // set order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHFFZVertex::doinit() {
  for(int ib=23;ib<34;ib+=10) {
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix,    ix,    ib);
    }
    addToList( -8,   8,  ib);
    addToList( -6,   8,  ib);
    addToList( -8,   6,  ib);
    // the leptons
    for(int ix=11;ix<17;++ix) {
      addToList(-ix,    ix,   ib);
    }
  }
  FFVVertex::doinit();
  cLHModelPtr model = dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHModel "
				   << " in LHFFZVertex::doinit()"
				   << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double pre =-0.5/sw/cw;
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double sp2(sqr(sp)),cp2(sqr(cp));
  // from Eqn A35
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  double xB(-2.5/sw*sp*cp*(cp2-sp2));
  double yu  = -0.4, ye  =  0.6;
  double vf(model->vev()/model->f());
  double xL(sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2())));
  double vu  = pre*( 0.5-4./3.*sw2-sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+7./15.  -cp2/6.)));
  double vd  = pre*(-0.5+2./3.*sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+11./15. +cp2/6.))); 
  double ve  = pre*(-0.5+2.*   sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*ye-9./5. +1.5*cp2))); 
  double vv  = pre*(+0.5          -sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(   ye-4./5. +0.5*cp2)));
  double au  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(0.2-0.5*cp2)));
  double ad  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s-sw*xB/sp/cp*(0.2-0.5*cp2)));
  double ae  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s-sw*xB/sp/cp*(0.2-0.5*cp2)));
  double av  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(0.2-0.5*cp2)));
  double vtl = pre*( 0.5-4./3.*sw2-sqr(vf)*(-0.5*sqr(xL)+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+9./5.-1.5*cp2
							  +(7./15.-2.*cp2/3.)*xL)));
  double atl = pre*(-0.5          -sqr(vf)*(+0.5*sqr(xL)-0.5*cw*xW*c/s
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
  _gl[6] = vtl - atl;
  _gr[6] = vtl + atl;
  _gl[7] = vtm - atm;
  _gr[7] = vtm + atm;
  _gl[8] = vth - ath;
  _gr[8] = vth + ath;
  // heavy Z
  double fact = 0.25*c/s/sw;
  vu  =  fact;
  vd  = -fact;
  ve  = -fact;
  vv  =  fact;
  au  = -fact;
  ad  =  fact;
  ae  =  fact;
  av  = -fact;
  vtl =  fact;
  atl = -fact;
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
  _glH[6] = vtl - atl;
  _grH[6] = vtl + atl;
  _glH[7] = vtm - atm;
  _grH[7] = vtm + atm;
  _glH[8] = vth - ath;
  _grH[8] = vth + ath;
}

void LHFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm = abs(a->id());
  int ianti = abs(b->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)|| iferm == 8);
  // Z0
  if(c->id()==ParticleID::Z0) {
    if(ianti==iferm) {
      left (_gl[iferm]);
      right(_gr[iferm]);
    }
    else {
      left (_gl[7]);
      right(_gr[7]);
    }
  }
  else {
    if(ianti==iferm) {
      left (_glH[iferm]);
      right(_grH[iferm]);
    }
    else {
      left (_glH[7]);
      right(_grH[7]);
    }
  }
}
