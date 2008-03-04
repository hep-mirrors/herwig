// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFPVertex class.
//

#include "LHTPFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "LHTPModel.h"

using namespace Herwig;

void LHTPFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge << _coupd << _coupu << _coupe << _coupnu 
     << _tmtpL << _tmtpR << _tmtL << _tmtR;
}

void LHTPFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >> _coupd >> _coupu >> _coupe >> _coupnu 
     >> _tmtpL >> _tmtpR >> _tmtL >> _tmtR;
}

ClassDescription<LHTPFFPVertex> 
LHTPFFPVertex::initLHTPFFPVertex;
// Definition of the static class description member.

void LHTPFFPVertex::Init() {

  static ClassDocumentation<LHTPFFPVertex> documentation
    ("The LHTPFFPVertex class implements the coupling"
     " of the charged fermions to the photon in the Little Higgs"
     " model with T-parity.");

}

LHTPFFPVertex::LHTPFFPVertex() : 
  _charge(37,0.0), _couplast(0.), _q2last(-1.*GeV2),
  _coupd(0.), _coupu(0.), _coupe(0.), _coupnu(0.),
  _tmtpL(0.), _tmtpR(0.), _tmtL(0.), _tmtR(0.) {
  // PDG codes for the particles
  vector<int> first,second;
  // interactions with the photon
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
  }
//   // extra top quark
//   first.push_back(-8);
//   second.push_back(8);
//   // the T-odd quarks
//   for(unsigned int ix=4000001;ix<4000007;++ix) {
//     first.push_back(-ix);
//     second.push_back(ix);
//   }
//   // the T-odd leptons
//   for(unsigned int ix=4000011;ix<4000017;ix+=2) {
//     first.push_back(-ix);
//     second.push_back(ix);
//   }
//   // extra top quark
//   first.push_back(-4000008);
//   second.push_back(4000008);
  vector<int> third(first.size(),22);
//   // interactions with A_H
//   // quark and T-odd quark
//   for(unsigned int ix=1;ix<7;++ix) {
//     first.push_back(-ix-4000000);
//     second.push_back(ix);
//     third.push_back(32);
//     first.push_back(-ix);
//     second.push_back(ix+4000000);
//     third.push_back(32);
//   }
//   // leptons and T-odd leptons
//   for(unsigned int ix=11;ix<17;ix+=2) {
//     first.push_back(-ix-4000000);
//     second.push_back(ix);
//     third.push_back(32);
//     first.push_back(-ix);
//     second.push_back(ix+4000000);
//     third.push_back(32);
//   }
//   // T+T-A_H
//   first .push_back(-4000008);
//   second.push_back(       8);
//   third .push_back(      32);
//   first .push_back( 4000008);
//   second.push_back(      -8);
//   third .push_back(      32);
//   // T-tA_H
//   first .push_back(-4000008);
//   second.push_back(       6);
//   third .push_back(      32);
//   first .push_back( 4000008);
//   second.push_back(      -6);
//   third .push_back(      32);
  setList(first,second,third);
}

void LHTPFFPVertex::doinit() throw(InitException) {
  // charges
  for(unsigned int ix=1;ix<16;++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) _charge[ix] = double(ptemp->iCharge())/3.;
  }
  for(unsigned int ix=4000001;ix<4000016;++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) _charge[ix-3999980] = double(ptemp->iCharge())/3.;
  }
  // couplings to A_H
  double sw = generator()->standardModel()->sin2ThetaW();
  double cw = sqrt(1.-sw);
  sw = sqrt(sw);
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFPVertex::doinit()"
				   << Exception::runerror;
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  // couplings of fermion T-odd fermion A_H
  _coupd  = -0.1*(cH/cw-5.*sH/sw);
  _coupu  = -0.1*(cH/cw+5.*sH/sw);
  _coupe  = -0.1*(cH/cw-5.*sH/sw);
  _coupnu = -0.1*(cH/cw+5.*sH/sw);
  // couplings of T+T- A_H
  _tmtpL = 0.4*cH/cw*model->cosBeta();
  _tmtpR = 0.4*cH/cw*model->cosAlpha();
  // couplings of T-t A_H
  _tmtL  = 0.4*cH/cw*model->sinBeta();
  _tmtR  = 0.4*cH/cw*model->sinAlpha();
  orderInGem(1);
  orderInGs(0);
  FFVVertex::doinit();
}

// coupling for FFP vertex
void LHTPFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id()),ibos(c->id());
  if(ibos==ParticleID::gamma) {
    if((iferm<20)) {
      setLeft(_charge[iferm]);
      setRight(_charge[iferm]);
    }
    else {
      iferm-=3999980;
      setLeft(_charge[iferm]);
      setRight(_charge[iferm]);
    }
  }
  else {
    int ianti = abs(b->id())%4000000;
    iferm = iferm%4000000;
    if(iferm==8) {
      if(ianti==8) {
	setLeft (_tmtpL);
	setRight(_tmtpR);
      }
      else {
	setLeft (_tmtL);
	setRight(_tmtR);
      }
    }
    if(iferm<=6) {
      if(iferm%2==0) setLeft(_coupu );
      else           setLeft(_coupd );
    }
    else {
      if(iferm%2==0) setLeft(_coupnu);
      else           setLeft(_coupe );
    }
    setRight(0.);
  }
}
