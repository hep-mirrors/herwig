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

IBPtr LHTPFFPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFPVertex::fullclone() const {
  return new_ptr(*this);
}

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
  orderInGem(1);
  orderInGs(0);
  // interactions with the photon
  // the quarks
  for(unsigned int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 22);
  }
  // the leptons
  for(unsigned int ix = 11; ix < 17; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-8,   8, 22);
  // the T-odd quarks
  for(long ix = 4000001;ix < 4000007; ++ix) {
    addToList(-ix,     ix, 22);
  }
  // the T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-4000008,  4000008, 22);
  // interactions with A_H
  // quark and T-odd quark
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix - 4000000,    ix,    32);
    addToList(-ix,    ix + 4000000,    32);
  }
  // leptons and T-odd leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix - 4000000,    ix,    32);
    addToList(-ix,    ix + 4000000,    32);
  }
  // T+T-A_H
  addToList(-4000008,         8,        32);
  addToList( 4000008,        -8,        32);
  // T-tA_H
  addToList(-4000008,         6,        32);
  addToList( 4000008,        -6,        32);
}

void LHTPFFPVertex::doinit() {
  // charges
  for(int ix = 1; ix < 16; ++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) _charge[ix] = double(ptemp->iCharge())/3.;
  }
  for(int ix = 4000001; ix < 4000016; ++ix) {
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
  FFVVertex::doinit();
}

// coupling for FFP vertex
void LHTPFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  long iferm=abs(a->id()),ibos(c->id());
  if(ibos == ParticleID::gamma) {
    if(iferm < 20) {
      left(_charge[iferm]);
      right(_charge[iferm]);
    }
    else {
      iferm-=3999980;
      left(_charge[iferm]);
      right(_charge[iferm]);
    }
  }
  else {
    long ianti = abs(b->id()) % 4000000;
    iferm = iferm % 4000000;
    if(iferm == 8) {
      if(ianti == 8) {
	left (_tmtpL);
	right(_tmtpR);
      }
      else {
	left (_tmtL);
	right(_tmtR);
      }
    }
    if(iferm <= 6) {
      if(iferm % 2 == 0) left(_coupu );
      else left(_coupd );
    }
    else {
      if(iferm %2 == 0) left(_coupnu);
      else left(_coupe );
    }
    right(0.);
  }
}
