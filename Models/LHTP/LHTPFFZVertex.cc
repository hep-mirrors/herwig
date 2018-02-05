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
  os << gl_ << gr_ << tl_ << tr_ << coupd_ << coupu_ << coupe_ << coupnu_ 
     << sL_ << cL_ << sR_ << cR_;
}

void LHTPFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> gl_ >> gr_ >> tl_ >> tr_ >> coupd_ >> coupu_ >> coupe_ >> coupnu_
     >> sL_ >> cL_ >> sR_ >> cR_;
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

LHTPFFZVertex::LHTPFFZVertex() 
  : gl_(37,0.0), gr_(37,0.0), tl_( 6,0.0), tr_( 6,0.0),
    coupd_(0.), coupu_(0.), coupe_(0.), coupnu_(0.),
    sL_(0.), cL_(1.), sR_(0.), cR_(1.),
    coupLast_(0.0), q2Last_(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHTPFFZVertex::doinit() {
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
  for(long ix = 11;ix<17;++ix) {
    addToList(-ix-4000000,    ix+4000000,    23);
  }
  // Z_H
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }  
  addToList(      -8,  +4000008,  33);
  addToList(       8,  -4000008,  33);
  addToList(      -6,  +4000008,  33);
  addToList(       6,  -4000008,  33);
  // the leptons
  for(int ix=11;ix<17;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }
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
    gl_[2*ix-1]  = fact*(model->vd()  + model->ad() );
    gl_[2*ix ]   = fact*(model->vu()  + model->au() );
    gl_[2*ix+9 ] = fact*(model->ve()  + model->ae() );
    gl_[2*ix+10] = fact*(model->vnu() + model->anu());
    gr_[2*ix-1]  = fact*(model->vd()  - model->ad() );
    gr_[2*ix ]   = fact*(model->vu()  - model->au() );
    gr_[2*ix+9 ] = fact*(model->ve()  - model->ae() );
    gr_[2*ix+10] = fact*(model->vnu() - model->anu());
    // T-odd fermions
    gl_[2*ix-1 +20] = fact*(model->vd()  + model->ad() );
    gl_[2*ix   +20] = fact*(model->vu()  + model->au() );
    gl_[2*ix+9 +20] = fact*(model->ve()  + model->ae() );
    gl_[2*ix+10+20] = fact*(model->vnu() + model->anu());
    gr_[2*ix-1 +20] = gl_[2*ix-1 +20];
    gr_[2*ix   +20] = gl_[2*ix   +20];
    gr_[2*ix+9 +20] = gl_[2*ix+9 +20];
    gr_[2*ix+10+20] = gl_[2*ix+10+20];
  }
  // couplngis to Z for extended top sector
  tl_[0] = (0.5*sqr(model->cosThetaL())-2./3.*sqr(sw))/cw/sw;
  tr_[0] = -2./3.*sw/cw;
  tl_[1] = (0.5*sqr(model->sinThetaL())-2./3.*sqr(sw))/cw/sw;
  tr_[1] = -2./3.*sw/cw;
  tl_[2] =  0.5/sw/cw*model->sinThetaL()*model->cosThetaL();
  tr_[2] = 0.;
  // couplings to the Z_H of T-odd fermions
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  sL_ = model->sinThetaL();
  cL_ = model->cosThetaL();
  sR_ = model->sinThetaR();
  cR_ = model->cosThetaR();
  coupd_  = 0.1*(sH/cw+5.*cH/sw);
  coupu_  = 0.1*(sH/cw-5.*cH/sw);
  coupe_  = 0.1*(sH/cw+5.*cH/sw);
  coupnu_ = 0.1*(sH/cw-5.*cH/sw);
  tl_[5] = -2./3.*sw/cw;
  tr_[5] = -2./3.*sw/cw;
  // couplings of T-odd top to the Z_H
  tl_[3] = 0.4*sH*sL_/cw;
  tr_[3] = 0.4*sH*sR_/cw;
  tl_[4] =-0.4*sH*cL_/cw;
  tr_[4] =-0.4*sH*cR_/cw;
  // base class initialisation
  FFVVertex::doinit();
}

void LHTPFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -electroMagneticCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  long iferm = abs(a->id());
  long ianti = abs(b->id());
  long ibos  = c->id();
  if(ibos == ParticleID::Z0) {
    if(iferm == 8 || iferm == 6) {
      if(iferm == 6 && ianti == 6) {
	left (tl_[0]);
	right(tr_[0]);
      }
      else if(iferm == 8 && ianti == 8) {
	left (tl_[1]);
	right(tr_[1]);
      }
      else {
	left (tl_[2]);
	right(tr_[2]);
      }
    }
    else if(iferm == 4000008) {
      left (tl_[5]);
      right(tr_[5]);
    }
    else if((iferm >= 1 && iferm <= 6)|| (iferm >= 11 && iferm <= 16)) {
      left (gl_[iferm]);
      right(gr_[iferm]);
    }
    else {
      iferm = (iferm % 4000000) + 20;
      left (gl_[iferm]);
      right(gr_[iferm]);
    }
  }
  else if(ibos == 33) {
    if(iferm>4000000) swap(iferm,ianti);
    assert(iferm<4000000&&ianti>4000000);
    if( iferm == 6 || iferm == 8 ) {
      if     (iferm==6&&ianti==4000006) {
	left (cL_*coupu_);
	right(0.);
      }
      else if(iferm==6&&ianti==4000008) {
	left (tl_[3]);
	right(tr_[3]);
      }
      else if(iferm==8&&ianti==4000006) {
	left (sL_*coupu_);
	right(0.);
      }
      else if(iferm==8&&ianti==4000008) {
	left ( tl_[4]);
	right( tr_[4]);
      }
      else
	assert(false);
    }
    else {
      right(0.);
      if(iferm <= 6) {
	if(iferm % 2 == 0) left( coupu_ );
	else               left( coupd_ );
      }
      else {
	if(iferm % 2 == 0) left( coupnu_ );
	else               left( coupe_  );
      }
    }
  }
  else
    assert(false);
}
