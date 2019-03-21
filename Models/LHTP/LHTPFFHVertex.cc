// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFHVertex class.
//

#include "LHTPFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHTPFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(cL_,1./GeV) << ounit(cR_,1./GeV) << model_;
}

void LHTPFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(cL_,1./GeV) >> iunit(cR_,1./GeV) >> model_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFHVertex,FFSVertex>
describeHerwigLHTPFFHVertex("Herwig::LHTPFFHVertex", "HwLHTPModel.so");

void LHTPFFHVertex::Init() {

  static ClassDocumentation<LHTPFFHVertex> documentation
    ("The LHTPFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model with T-parity");

}

LHTPFFHVertex::LHTPFFHVertex() 
  : q2Last_(ZERO) {
  orderInGem(1);
  orderInGs(0);
  massLast_[0] = 0.*GeV; 
  massLast_[1] = 0.*GeV;
  idLast_[0] = 0;
  idLast_[1] = 0;
}

void LHTPFFHVertex::doinit() {
  // SM like higgs
  addToList(  -3,   3,  25);
  addToList(  -4,   4,  25);
  addToList(  -5,   5,  25);
  addToList(  -6,   6,  25);
  addToList(  -6,   8,  25);
  addToList(  -8,   6,  25);
  addToList(  -8,   8,  25);
  addToList( -13,  13,  25);
  addToList( -15,  15,  25);
  addToList( -4000002, 4000002,  25);
  addToList( -4000004, 4000004,  25);
  addToList( -4000006, 4000006,  25);
  addToList( -4000012, 4000012,  25);
  addToList( -4000014, 4000014,  25);
  addToList( -4000016, 4000016,  25);
  // phi0
  addToList( -3      , 4000003,  35);
  addToList( -4      , 4000004,  35);
  addToList( -5      , 4000005,  35);
  addToList( -4000003,       3,  35);
  addToList( -4000004,       4,  35);
  addToList( -4000005,       5,  35);
  addToList( -6      , 4000006,  35);
  addToList( -8      , 4000006,  35);
  addToList( -4000006,       6,  35);
  addToList( -4000006,       8,  35);
  // phiP
  addToList( -2      , 4000002,  36);
  addToList( -3      , 4000003,  36);
  addToList( -4      , 4000004,  36);
  addToList( -5      , 4000005,  36);
  addToList( -4000002,       2,  36);
  addToList( -4000003,       3,  36);
  addToList( -4000004,       4,  36);
  addToList( -4000005,       5,  36);
  addToList( -12     , 4000012,  36);
  addToList( -14     , 4000014,  36);
  addToList( -16     , 4000016,  36);
  addToList( -4000012,      12,  36);
  addToList( -4000014,      14,  36);
  addToList( -4000016,      16,  36);
  addToList( -6      , 4000006,  36);
  addToList( -6      , 4000008,  36);
  addToList( -8      , 4000006,  36);
  addToList( -8      , 4000008,  36);
  addToList( -4000006,       6,  36);
  addToList( -4000008,       6,  36);
  addToList( -4000006,       8,  36);
  addToList( -4000008,       8,  36);
  // phi +/-
  addToList( -1      , 4000002, -37);
  addToList( -3      , 4000004, -37);
  addToList( -5      , 4000006, -37);
  addToList( -4000001,       2, -37);
  addToList( -4000003,       4, -37);
  addToList( -4000005,       6, -37);
  addToList( -4000005,       8, -37);
  addToList( -4000002,       1,  37);
  addToList( -4000004,       3,  37);
  addToList( -4000006,       5,  37);
  addToList( -2      , 4000001,  37);
  addToList( -4      , 4000003,  37);
  addToList( -6      , 4000005,  37);
  addToList( -8      , 4000005,  37);
  addToList( -11     , 4000012, -37);
  addToList( -13     , 4000014, -37);
  addToList( -15     , 4000016, -37);
  addToList( -4000011,      12, -37);
  addToList( -4000013,      14, -37);
  addToList( -4000015,      16, -37);
  addToList( -4000012, 11     ,  37);
  addToList( -4000014, 13     ,  37);
  addToList( -4000016, 15     ,  37);
  addToList(      -12, 4000011,  37);
  addToList(      -14, 4000013,  37);
  addToList(      -16, 4000015,  37);
  model_ = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model_)   throw InitException() << "Must be using the LHModel "
				      << " in LHFFPVertex::doinit()"
				      << Exception::runerror;
  cL_  .resize(18); 
  cR_  .resize(18);
  Energy v  = model_->vev();
  Energy f  = model_->f();
  double vf = model_->vev()/model_->f();
  double sa = model_->sinAlpha();
  double ca = model_->cosAlpha();
  // lightest higgs couplings
  // coupling of light SM fermions
  cL_[0]   =  cR_[0] = 1./v;
  // couplings to top quarks
  cL_[1]   =  cR_[1] = sa*ca/f;
  cL_[2]   = -sa/ca/v;
  cR_[2]   =  sqr(ca)*vf/v;
  // couplings to T-odd quarks
  cL_[3]   =  cR_[3] = 0.5*sqrt(0.5)/f*model_->kappaQuark();
  // couplings to T-odd leptons
  cL_[4]   =  cR_[4] = 0.5*sqrt(0.5)/f*model_->kappaLepton();
  // Phi0
  // quark, T-odd quark
  cL_[5]   = sqrt(0.5)/f;
  cR_[5]   = ZERO;
  // and top quarks
  cL_[6]   = sqrt(0.5)*model_->cosThetaR()/f/ca;
  cR_[6]   = ZERO;
  cL_[7]   = sqrt(0.5)*model_->sinThetaR()/f/ca;
  cR_[7]   = ZERO;
  // PhiP
  // quark, T-odd quark
  cL_[8] = vf/f*model_->kappaQuark()/12;
  cR_[8] = sqrt(0.5)/f;
  // lepton, T-odd lepton
  cL_[9] = vf/f*model_->kappaLepton()/12;
  cR_[9] = ZERO;
  // top, T_-
  cL_[10] = model_->cosThetaR()*sqrt(2.)*vf/f/ca/3;
  cR_[10] = ZERO;
  cL_[11] = model_->sinThetaR()*sqrt(2.)*vf/f/ca/3;
  cR_[11] = ZERO;
  // top, t_-
  cL_[12] = model_->cosThetaL()*vf/f/12.*model_->kappaQuark();
  cR_[12] = model_->cosThetaR()*sqrt(0.5)/f/ca;
  cL_[13] = model_->sinThetaL()*vf/f/12.*model_->kappaQuark();
  cR_[13] = model_->sinThetaR()*sqrt(0.5)/f/ca;
  // Phi +/-
  cL_[14] = vf/f*model_->kappaLepton()/24.;
  cR_[14] = ZERO;
  // quark T-odd quark
  cL_[15] = vf/f*model_->kappaQuark() /24.;
  cR_[15] =-vf*sqrt(0.5)/v;
  cL_[16] = model_->cosThetaL()*vf/f*model_->kappaQuark() /24.;
  cR_[16] =-model_->cosThetaR()*vf*sqrt(0.5)/v/ca;
  cL_[17] = model_->sinThetaL()*vf/f*model_->kappaQuark() /24.;
  cR_[17] =-model_->sinThetaR()*vf*sqrt(0.5)/v/ca;
  FFSVertex::doinit();
}

void LHTPFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  norm(1.);
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  //  int ihigg=abs(c->id());
  // left and right couplings set to one
  // SM like higgs
  if(c->id()==ParticleID::h0) {
    // to SM fermions and T
    if(iferm<=16&&ianti<=16) {
      // running masses
      if(q2!=q2Last_||idLast_[0]!=iferm||idLast_[1]!=ianti) {
	q2Last_ = q2;
	idLast_[0] = iferm;
	assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	       (idLast_[0]>=11 && idLast_[0]<=16));
	if(idLast_[0]!=8)
	  massLast_[0] = model_->mass(q2,a);
	else
	  massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
	idLast_[1] = ianti;
	assert((idLast_[1]>=1  && idLast_[1]<=8 ) || 
	       (idLast_[1]>=11 && idLast_[1]<=16));
	if(idLast_[0]!=idLast_[1]) {
	  if(idLast_[1]!=8)
	    massLast_[1] = model_->mass(q2,a);
	  else
	    massLast_[1] = model_->mass(q2,getParticleData(ParticleID::t));
	}
	else {
	  massLast_[1] = massLast_[0];
	}
      }
      if(iferm<6||iferm>8) {
	left (-Complex(cL_[0]*massLast_[0]));
	right(-Complex(cR_[0]*massLast_[0]));
      }
      else {
	if(iferm==8&&ianti==8) {
	  left ( Complex(cL_[1]*massLast_[0]));
	  right( Complex(cR_[1]*massLast_[0]));
	}
	else {
	  if(a->id()==ParticleID::tbar||b->id()==ParticleID::tbar) {
	    left (-Complex(cL_[2]*massLast_[0]));
	    right(-Complex(cR_[2]*massLast_[0]));
	  }
	  else {
	    left (-Complex(cR_[2]*massLast_[0]));
	    right(-Complex(cL_[2]*massLast_[0]));
	  } 
	}
      }
    }
    else {
      if(iferm<=4000006) {
	left ( Complex(cL_[3]*model_->vev()));
	right( Complex(cR_[3]*model_->vev()));
      }
      else {
	left ( Complex(cL_[4]*model_->vev()));
	right( Complex(cR_[4]*model_->vev()));
      }
    }
  }
  // Phi0
  else if(c->id()==ParticleID::H0 ||
	  c->id()==ParticleID::A0) {
    tcPDPtr ferm = a;
    if(iferm>4000000) {
      swap(iferm,ianti);
      ferm = b;
    }
    if(q2!=q2Last_||idLast_[0]!=iferm) {
      q2Last_ = q2;
      idLast_[0] = iferm;
      assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	     (idLast_[0]>=11 && idLast_[0]<=16));
      if(idLast_[0]!=8)
	massLast_[0] = model_->mass(q2,ferm);
      else
	massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
    }
    if(c->id()==ParticleID::H0 ) {
      unsigned int      iloc = 5;
      if(iferm==6)      iloc = 6;
      else if(iferm==8) iloc = 7;
      if( (a->id()>=1&&a->id()<=8) || (b->id()>=1&&b->id()<=8) ) {
	left ( Complex(cR_[iloc]*massLast_[0]));
	right( Complex(cL_[iloc]*massLast_[0]));
      }
      else {
	left ( Complex(cL_[iloc]*massLast_[0]));
	right( Complex(cR_[iloc]*massLast_[0]));
      }
    }
    // PhiP
    else if(c->id()==ParticleID::A0) {
      if(iferm<=5) {
	if( (a->id()>=1&&a->id()<=5) || (b->id()>=1&&b->id()<=5) ) {
	  if(iferm%2==0) {
	    right(Complex(0., 1.)*model_->vev()*cL_[8]);
	    left (Complex(0.,-1.)*massLast_[0] *cR_[8]);
	  }
	  else {
	    right(Complex(ZERO));
	    left (Complex(0., 1.)*massLast_[0] *cR_[8]);
	  }
	}
	else {
	  if(iferm%2==0) {
	    left (Complex(0.,-1.)*model_->vev()*cL_[8]);
	    right(Complex(0., 1.)*massLast_[0] *cR_[8]);
	  }
	  else {
	    left (Complex(ZERO));
	    right(Complex(0.,-1.)*massLast_[0] *cR_[8]);
	  }
	}
      }
      else if(iferm>=12) {
	if( (a->id()>=11&&a->id()<=16) || (b->id()>=11&&b->id()<=16) ) {
	  right(Complex(0., 1.)*model_->vev()*cL_[9]);
	  left (Complex(0.,-1.)*massLast_[0] *cR_[9]);
	}
	else {
	  right(Complex(0.,-1.)*massLast_[0] *cR_[9]);
	  left (Complex(0., 1.)*model_->vev()*cL_[9]);
	}
      }
      else {
	if(ianti==4000008) {
	  unsigned int iloc = (iferm+14)/2;
	  if( (a->id()==6||a->id()==8) || (b->id()==6||b->id()==8) ) {
	    left (Complex(0., 1.)*massLast_[0]*cR_[iloc]);
	    right(Complex(0., 1.)*massLast_[0]*cL_[iloc]);
	  }
	  else {
	    left (Complex(0.,-1.)*massLast_[0]*cL_[iloc]);
	    right(Complex(0.,-1.)*massLast_[0]*cR_[iloc]);
	  }
	}
	else {
	  unsigned int iloc = (iferm+18)/2;
	  if( (a->id()==6||a->id()==8) || (b->id()==6||b->id()==8) ) {
	    left (Complex(0., 1.)*model_->vev()*cL_[iloc]);
	    right(Complex(0.,-1.)*massLast_[0] *cR_[iloc]);
	  }
	  else {
	    left (Complex(0., 1.)*massLast_[0] *cR_[iloc]);
	    right(Complex(0.,-1.)*model_->vev()*cL_[iloc]);
	  }
	}
      }
    }
  }
  else if(abs(c->id())==ParticleID::Hplus) {
    tcPDPtr ferm = a;
    if(iferm>4000000) {
      swap(iferm,ianti);
      ferm = b;
    }
    if(q2!=q2Last_||idLast_[0]!=iferm) {
      q2Last_ = q2;
      idLast_[0] = iferm;
      assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	     (idLast_[0]>=11 && idLast_[0]<=16));
      if(idLast_[0]!=8)
	massLast_[0] = model_->mass(q2,ferm);
      else
	massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
    }
    Complex cleft(0.),cright(0.);
    // lepton and T-odd lepton
    if(iferm>=11&&iferm<=16) {
      cright = Complex(cR_[14]*massLast_[0]);
      cleft  = Complex(cL_[14]*model_->vev()); 
    }
    else if(iferm>=1&&iferm<=6) {
      cright = Complex(cR_[15]*massLast_[0]);
      cleft  = Complex(cL_[15]*model_->vev());
    }
    else if(iferm==6) {
      cright = Complex(cR_[16]*massLast_[0]);
      cleft  = Complex(cL_[16]*model_->vev()); 
    }
    else if(iferm==8) {
      cright = Complex(cR_[17]*massLast_[0]);
      cleft  = Complex(cL_[17]*model_->vev()); 
    }
    if((a->id()>=1&&a->id()<=16) ||(b->id()>=1&&b->id()<=16) ) {
      swap(cleft,cright);
      cleft  *= -1.;
      cright *= -1.;
    }
    if(c->id()==ParticleID::Hminus) {
      cleft  *= -1.;
      cright *= -1.;
    }
    left (cleft );
    right(cright);
  }
}
