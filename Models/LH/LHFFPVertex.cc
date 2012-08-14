// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFPVertex class.
//

#include "LHFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHFFPVertex::LHFFPVertex() 
  : _couplast(0.), _q2last(-1.*GeV2) {
  // order in strong and em coupling
  orderInGem(1);
  orderInGs(0);
}

void LHFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge <<  _gl << _gr;
}

void LHFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >>  _gl >> _gr;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFPVertex,FFVVertex>
describeHerwigLHFFPVertex("Herwig::LHFFPVertex", "HwLHModel.so");

void LHFFPVertex::Init() {

  static ClassDocumentation<LHFFPVertex> documentation
    ("The LHFFPVertex class implements the couplings of the fermions to"
     " the photon and A_H in the Little Higgs model");

}

void LHFFPVertex::doinit() {
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,    ix,    22);
    addToList(-ix,    ix,    32);
  }
  addToList( -8,   8,  22);
  addToList( -8,   8,  32);
  addToList( -6,   8,  32);
  addToList( -8,   6,  32);
  // the leptons
  for(int ix=11;ix<17;++ix) {
    if(ix%2!=0) addToList(-ix,    ix,    22);
    addToList(-ix,    ix,    32);
  }
  FFVVertex::doinit();
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHFFPVertex::doinit()"
			  << Exception::runerror;
  // charges
  _charge.resize(17);  
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = model->ed();
    _charge[2*ix ]   = model->eu();
    _charge[2*ix+9 ] = model->ee();
    _charge[2*ix+10] = model->enu();
  }
  _charge[8] =  model->eu();
  // couplings for the heavy photon taken from table IX
  double cw  = sqrt(1.-sin2ThetaW());
  double xL  = sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2()));
  double cp2 = sqr(model->cosThetaPrime());
  double yu  = -0.4;
  double ye  =  0.6;
  // prefactor after removal of -e
  double pre = -0.5/cw/model->cosThetaPrime()/model->sinThetaPrime();
  // down type quarks
  double gvd   = pre*(2.*yu+11./15.+cp2/6);
  double gad   = pre*(-0.2+0.5*cp2);
  // up type quarks
  double gvu   = pre*(2.*yu+17./15.-5./6.*cp2);
  double gau   = pre*( 0.2-0.5*cp2);
  // charged leptons
  double gve   = pre*(2.*ye-9./5.+1.5*cp2);
  double gae   = pre*(-0.2+0.5*cp2);
  // neutrinos
  double gvv   = pre*(-0.2+0.5*cp2);
  double gav   = pre*( 0.2-0.5*cp2);
  // light top
  double gvtll = pre*(2.*yu+17./15.-5./6.*cp2-0.2*xL);
  double gatll = pre*(0.2-0.5*cp2-0.2*xL);
  // mixed top
  double gvtlh = pre*0.2*model->lambda1()*model->lambda2()/
    (sqr(model->lambda1())+sqr(model->lambda2()));
  double gatlh = gvtlh;
  // heavy top
  double gvthh = pre*(2.*yu+14./15.-4./3.*cp2+0.2*xL);
  double gathh = pre*0.2*xL;
  _gl.resize(17);
  _gr.resize(17);
  for(unsigned int ix=1;ix<4;++ix) {
    _gr[2*ix-1]  = gvd+gad;
    _gl[2*ix-1]  = gvd-gad;
    _gr[2*ix ]   = gvu+gau;
    _gl[2*ix ]   = gvu-gau;
    _gr[2*ix+9 ] = gve+gae;
    _gl[2*ix+9 ] = gve-gae;
    _gr[2*ix+10] = gvv+gav;
    _gl[2*ix+10] = gvv-gav;
  }
  // light top
  _gr[6] = gvtll+gatll;
  _gl[6] = gvtll-gatll;
  // mixed top
  _gr[7] = gvtlh+gatlh;
  _gl[7] = gvtlh-gatlh;
  // heavy top
  _gr[8] = gvthh+gathh;
  _gl[8] = gvthh-gathh;
}

// coupling for FFP vertex
void LHFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int iferm=abs(a->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)||iferm==8);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  // photon
  if(c->id()==ParticleID::gamma) {
    left (_charge[iferm]);
    right(_charge[iferm]);
  }
  // heavy photon
  else {
    assert(c->id()==32);
    int ianti = abs(b->id());
    if(ianti==iferm) {
      left (_gl[iferm]);
      right(_gr[iferm]);
    }
    else {
      left (_gl[7]);
      right(_gr[7]);
    }
  }
}
