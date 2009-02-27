// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFHVertex class.
//

#include "LHFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coup,1./GeV) << _theSM;
}

void LHFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,1./GeV) >> _theSM;
}

ClassDescription<LHFFHVertex> LHFFHVertex::initLHFFHVertex;
// Definition of the static class description member.

void LHFFHVertex::Init() {

  static ClassDocumentation<LHFFHVertex> documentation
    ("The LHFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model");

}

LHFFHVertex::LHFFHVertex() 
  : _idlast(0), _q2last(0.*GeV2) {
  _masslast[0] = 0.*GeV; 
  _masslast[1] = 0.*GeV; 
  // PDG codes for the particles
  vector<long> first,second,third;
  // the quarks
  for(int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  // the quarks
  for(int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(35);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(35);
  }
  // the quarks
  for(int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(36);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(36);
  }
  setList(first,second,third);
}

void LHFFHVertex::doinit() {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!_theSM) throw InitException() << "Must be using the LHModel "
				     << " in LHFFPVertex::doinit()"
				     << Exception::runerror;;
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model)   throw InitException() << "Must be using the LHModel "
				     << " in LHFFPVertex::doinit()"
				     << Exception::runerror;
  _coup.resize(10);
  Energy v   = model->vev();
  double s0  = model->sinTheta0();
  double s02 = sqr(s0);
  double vf  = model->vev()/model->f();
  double xL  = sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2()));
  double xR  = sqr(model->lambda1())/sqrt(sqr(model->lambda1())+sqr(model->lambda2()));
  Energy mT  = getParticleData(8)->mass();
  // coupling of light SM fermions
  _coup[0] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf))/v;
  // couplings to top quark
  _coup[1] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf)+sqr(vf)*xL*(1.+xL))/v;
  _coup[2] = xR*(1.+xL)*vf/mT;



  orderInGem(1);
  orderInGs(0);
  FFSVertex::doinit();
}

void LHFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c,
				       int) {
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  int ihigg=abs(c->id());
  // left and right couplings set to one
  setLeft(1.); setRight(1.);
  // first the overall normalisation
  if(q2!=_q2last) {
    _q2last=q2;
    _idlast=iferm;
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      _masslast[0]=_theSM->mass(q2,a);
    }
    if(iferm!=ianti) {
      if((ianti>=1 && ianti<=6)||(ianti>=11 &&ianti<=16)) {
	_masslast[1]=_theSM->mass(q2,b);
      }
    }
  }

//     else {
//       throw HelicityConsistencyError() << "LHFFHVertex::setCoupling " 
// 				       << "Unknown particle in Higgs vertex" 
// 				       << Exception::warning;
//       _masslast = 0*MeV;
//     }
//   }
//   else if(iferm!=_idlast) {
//     _idlast=iferm;
//     if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
//       _masslast=_theSM->mass(q2,a);
//     }
//     else {
//       throw HelicityConsistencyError() << "LHFFHVertex::setCoupling " 
// 				       << "Unknown particle in Higgs vertex" 
// 				       << Exception::warning;
//       _masslast = 0*MeV;
//     }
//   }
//   setNorm(_couplast*_masslast);
}
