// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsFFPVertex class.
//

#include "LittleHiggsFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

inline LittleHiggsFFPVertex::LittleHiggsFFPVertex() 
  : _couplast(0.), _q2last(-1.*GeV2) {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(22);
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(32);
  }
  first.push_back( -8);
  second.push_back( 8);
  third.push_back(22);
  first.push_back( -8);
  second.push_back( 8);
  third.push_back(32);
  first.push_back( -6);
  second.push_back( 8);
  third.push_back(32);
  first.push_back( -8);
  second.push_back( 6);
  third.push_back(32);
  // the leptons
  for(unsigned int ix=11;ix<17;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(22);
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(32);
  }
  setList(first,second,third);
}

void LittleHiggsFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge <<  _gv << _ga << _model;
}

void LittleHiggsFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >>  _gv >> _ga >> _model;
}

ClassDescription<LittleHiggsFFPVertex> LittleHiggsFFPVertex::initLittleHiggsFFPVertex;
// Definition of the static class description member.

void LittleHiggsFFPVertex::Init() {

  static ClassDocumentation<LittleHiggsFFPVertex> documentation
    ("There is no documentation for the LittleHiggsFFPVertex class");

}

void LittleHiggsFFPVertex::doinit() throw(InitException) {
  FFVVertex::doinit();
  _model = dynamic_ptr_cast<cLittleHiggsModelPtr>(generator()->standardModel());
  if(!_model) 
    throw InitException() << "Must be using the LittleHiggsModel "
			  << " in LittleHiggsFFPVertex::doinit()"
			  << Exception::runerror;
  // charges
  _charge.resize(17);  
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = _model->ed();
    _charge[2*ix ]   = _model->eu();
    _charge[2*ix+9 ] = _model->ee();
    _charge[2*ix+10] = _model->enu();
  }
  _charge[8] =  _model->eu();
  // couplings for the heavy photon
  double cw  = sqrt(1.-_model->sin2ThetaW());
  double pre = 0.5/cw/_model->cosThetaPrime()/_model->sinThetaPrime();
  double xL  = sqr(_model->lambda1())/(sqr(_model->lambda1())+sqr(_model->lambda2()));
  double cp2 = sqr(_model->cosThetaPrime());
  double yu  = -0.4;
  double ye  =  0.6;
  // down type quarks
  double gvd   = pre*(2.*yu+11./15.+1./6.*cp2);
  double gad   = pre*(-0.2+0.5*cp2);
  // up type quarks
  double gvu   = pre*(2.*yu+17./15.-5./6.*cp2);
  double gau   = pre*( 0.2-0.5*cp2);
  // charged leptons
  double gve   = pre*(2.*ye-9./5.+1.5*cp2);
  double gae   = pre*(-0.2+0.5*cp2);
  // neutrinos
  double gvv   = pre*(ye-0.8+0.5*cp2);
  double gav   = pre*( 0.2-0.5*cp2);
  // light top
  double gvtll = pre*(2.*yu+17./15.-5./6.*cp2-0.2*xL);
  double gatll = pre*(0.2-0.5*cp2-0.2*xL);
  // mixed top
  double gvtlh = pre*0.2*_model->lambda1()*_model->lambda2()/
    (sqr(_model->lambda1())+sqr(_model->lambda2()));
  double gatlh = gvtlh;
  // heavy top
  double gvthh = pre*(2.*yu+14./15.-4./3.*cp2+0.2*xL);
  double gathh = pre*0.2*xL;
  _ga.resize(17);
  _gv.resize(17);
  for(unsigned int ix=1;ix<4;++ix) {
    _gv[2*ix-1]  = gvd;
    _ga[2*ix-1]  = gad;
    _gv[2*ix ]   = gvu;
    _ga[2*ix ]   = gau;
    _gv[2*ix+9 ] = gve;
    _ga[2*ix+9 ] = gae;
    _gv[2*ix+10] = gvv;
    _ga[2*ix+10] = gav;
  }
  // light top
  _gv[6] = gvtll;
  _ga[6] = gatll;
  // mixed top
  _gv[7] = gvtlh;
  _ga[7] = gatlh;
  // heavy top
  _gv[8] = gvthh;
  _ga[8] = gathh;
  // order in strong and em coupling
  orderInGem(1);
  orderInGs(0);
}

// coupling for FFP vertex
void LittleHiggsFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _model->alphaEM(q2);
    _couplast = -sqrt(4.0*Constants::pi*alpha);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)||iferm==8) {
    // photon
    if(c->id()==ParticleID::gamma) {
      setLeft(_charge[iferm]);
      setRight(_charge[iferm]);
    }
    // heavy photon
    else {
      int ianti = abs(b->id());
      if(ianti==iferm) {
	setRight(_gv[iferm]+_ga[iferm]);
	setLeft (_gv[iferm]-_ga[iferm]);
      }
      else {
	setRight(_gv[7]+_ga[7]);
	setLeft (_gv[7]-_ga[7]);
      }
    }
  }
  else
    throw HelicityConsistencyError() << "SMGFFPVertex::setCoupling "
				     << "Unknown particle in photon vertex" 
				     << Exception::runerror;
}
