// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWWWVertex class.
//

#include "LHWWWWVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWWWWVertex::LHWWWWVertex() : 
  _couplast(0.0), _q2last(sqr(Constants::MaxEnergy)), _coup(36,0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs(0);
}

IBPtr LHWWWWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWWWWVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _coup;
}

void LHWWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _coup;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWWWVertex,VVVVVertex>
describeHerwigLHWWWWVertex("Herwig::LHWWWWVertex", "HwLHModel.so");

void LHWWWWVertex::Init() {

  static ClassDocumentation<LHWWWWVertex> documentation
    ("The LHWWWWVertex class implements the quartic electroweak"
     " boson couplings in the Little Higgs Model");

}

void LHWWWWVertex::doinit() {
  // all charge W's
  addToList(24, -24, 24, -24);
  addToList(34, -34, 34, -34);
  addToList(24, -24, 34, -34);
  addToList(24, -24, 24, -34);
  addToList(24, -24, 34, -24);
  addToList(34, -24, 34, -24);
  addToList(24, -34, 24, -34);
  addToList(34, -34, 24, -34);
  addToList(34, -34, 34, -24);
  // two neutral and 2 W_L
  addToList(22,  24, 22, -24);
  addToList(23,  24, 23, -24);
  addToList(22,  24, 23, -24);
  addToList(22,  24, 32, -24);
  addToList(22,  24, 33, -24);
  addToList(23,  24, 33, -24);
  addToList(23,  24, 32, -24);
  addToList(33,  24, 33, -24);
  addToList(33,  24, 32, -24);
  // two neutral and 2 W_H
  addToList(22,  34, 22, -34);
  addToList(23,  34, 23, -34);
  addToList(22,  34, 23, -34);
  addToList(22,  34, 32, -34);
  addToList(22,  34, 33, -34);
  addToList(23,  34, 33, -34);
  addToList(23,  34, 32, -34);
  addToList(33,  34, 33, -34);
  addToList(33,  34, 32, -34);
  // two neutral W_L W_H
  addToList(23,  24, 23, -34);
  addToList(23,  24, 22, -34);
  addToList(22,  24, 32, -34);
  addToList(23,  24, 32, -34);
  addToList(33,  24, 33, -34);
  addToList(33,  24, 32, -34);
  addToList(22,  24, 33, -34);
  addToList(23,  24, 33, -34);
  addToList(23,  34, 23, -24);
  addToList(23,  34, 22, -24);
  addToList(22,  34, 32, -24);
  addToList(23,  34, 32, -24);
  addToList(33,  34, 33, -24);
  addToList(33,  34, 32, -24);
  addToList(22,  34, 33, -24);
  addToList(23,  34, 33, -24);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  // correction factors for the different interactions
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double vf(sqr(model->vev()/model->f()));
  double xB(-2.5/sw*sp*cp*(sqr(cp)-sqr(sp)));
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  double xH(2.5/sw/cw*s*c*sp*cp*(sqr(c*sp)+sqr(s*cp))/
	    (5.*sqr(sp*cp/sw)-sqr(s*c/cw)));
  // 4 W's
  _coup[ 0] =-1./sw2;
  _coup[ 1] =-1./sw2;
  _coup[ 2] = 0.5/sw2*(sqr(c)-sqr(s))/c/s;
  _coup[ 3] =-0.25/sw2*vf*s*c*(sqr(c)-sqr(s));
  _coup[ 4] =-0.25/sw2;
  _coup[ 5] =-1./sw2*(pow(c,6)+pow(s,c))/sqr(s*c);
  // 2 W_L
  _coup[ 6] = 1.;
  _coup[ 7] = sqr(cw/sw);
  _coup[ 8] = cw/sw;
  _coup[ 9] =-cw/sw*xB*vf;
  _coup[10] =-cw/sw*xW*vf+0.5/sw*s*c*(sqr(c)-sqr(s))*vf;
  _coup[11] =-(sqr(cw)-sw2)/sw2*xW*vf;
  _coup[12] =-sqr(cw/sw)*xB*vf;
  _coup[13] = 0.;
  _coup[14] = 1./sw2;
  _coup[15] = xH*vf/sw2;
  _coup[16] = 1.;
  _coup[17] = sqr(cw/sw);
  _coup[18] = cw/sw;
  _coup[19] =-cw/sw*xB*vf-xH/sw*vf*(sqr(c)-sqr(s))/s/c;
  _coup[20] =-1./sw*(sqr(c)-sqr(s))/s/c;
  _coup[21] =-cw/sw2*(sqr(c)-sqr(s))/s/c;
  _coup[22] =-sqr(cw/sw)*xB*vf-cw/sw2*xH*vf*(sqr(c)-sqr(s))/c/s;
  _coup[23] = 0.;
  _coup[24] = (pow(c,6)+pow(s,6))/sqr(s*c)/sw2;
  _coup[25] = xH/sw2*vf*(pow(c,6)+pow(s,6))/sqr(s*c)
             +cw/sw2*xB*vf*(sqr(c)-sqr(s))/s/c;
  _coup[26] = 0.;
  _coup[27] = 2.*cw/sw2*xW*vf;
  _coup[28] = 0.;
  _coup[29] = xH*vf/sw;
  _coup[30] = xH*vf*cw/sw2;
  _coup[31] = xW*vf/sw;
  _coup[32] =-(sqr(c)-sqr(s))/s/c/sw2;
  _coup[33] =-xH*vf*(sqr(c)-sqr(s))/s/c-cw/sw2*xB*vf;
  _coup[34] = 1./sw;
  _coup[35] = cw/sw2;
}

void LHWWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
			       tcPDPtr c,tcPDPtr d) {
  // id's of the particles
  long int id[4]={a->id(),b->id(),c->id(),d->id()};
  // order the particles
  int ngamma(0),nz(0);
  int iorder[4];
  for(int ix=0;ix<4;++ix) {
    if      (id[ix]==22||id[ix]==32) ++ngamma;
    else if (id[ix]==23||id[ix]==33) ++nz;
  }
  // if photons or Z's
  if(ngamma!=0 || nz!=0) {
    int iy=0;
    // put the photons first
    for(int ix=0;iy<ngamma&&ix<4;++ix) {
      if(id[ix]==22||id[ix]==32) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the Z bosons
    for(int ix=0;iy<ngamma+nz&&ix<4;++ix) {
      if(id[ix]==23||id[ix]==33) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24||id[ix]==34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==3);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24||id[ix]==-34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==4);
  }
  else {
    int iy=0;
    // first the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24||id[ix]==34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==2);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24||id[ix]==-34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==4);
  }
  setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);
  setType(2);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = sqr(electroMagneticCoupling(q2));
    _q2last=q2;
  }
  // ids of the particles
  for(unsigned int ix=0;ix<4;++ix) {
    if     (iorder[ix]==0) id[ix] = abs(a->id());
    else if(iorder[ix]==1) id[ix] = abs(b->id());
    else if(iorder[ix]==2) id[ix] = abs(c->id());
    else if(iorder[ix]==3) id[ix] = abs(d->id());
  }
  if( ngamma == 0 && nz == 0 ) {
    if(id[0]==id[1]) {
      if(id[2]==id[3]) {
	if(id[0]==24&&id[2]==24)
	  norm(_couplast*_coup[0]);
	else if(id[0]==34&&id[2]==34)
	  norm(_couplast*_coup[5]);
	else
	  norm(_couplast*_coup[1]);
      }
      else {
	if(id[0]==24)
	  norm(_couplast*_coup[3]);
	else
	  norm(_couplast*_coup[2]);
      }
    }
    else {
      if(id[2]==id[3]) {
	if(id[2]==24)
	  norm(_couplast*_coup[3]);
	else
	  norm(_couplast*_coup[2]);
      }
      else
	norm(_couplast*_coup[4]);
    } 
  }
  else {
    if(id[2]==id[3]) {
      unsigned int ioff = id[2]==24 ? 0 : 10;
      if(id[0]==22&&id[1]==22)
	norm(_couplast*_coup[6+ioff]);
      else if(id[0]==23&&id[1]==23)
	norm(_couplast*_coup[7+ioff]);
      else if((id[0]==22&&id[1]==23) || (id[0]==23&&id[1]==22))
	norm(_couplast*_coup[8+ioff]);
      else if((id[0]==22&&id[1]==32) || (id[0]==32&&id[1]==22))
	norm(_couplast*_coup[9+ioff]);
      else if((id[0]==22&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[10+ioff]);
      else if((id[0]==23&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[11+ioff]);
      else if((id[0]==23&&id[1]==32) || (id[0]==32&&id[1]==23))
	norm(_couplast*_coup[12+ioff]);
      else if( id[0]==33&&id[1]==33)
	norm(_couplast*_coup[14+ioff]);
      else if((id[0]==32&&id[1]==33) || (id[0]==33&&id[1]==32))
	norm(_couplast*_coup[15+ioff]);
      else
	assert(false);
    }
    else { 
      if(id[0]==23&&id[1]==23)
	norm(_couplast*_coup[27]);
      else if((id[0]==22&&id[1]==23) || (id[0]==23&&id[1]==22))
	norm(_couplast*_coup[28]);
      else if((id[0]==22&&id[1]==32) || (id[0]==32&&id[1]==22))
	norm(_couplast*_coup[29]);
      else if((id[0]==23&&id[1]==32) || (id[0]==32&&id[1]==23))
	norm(_couplast*_coup[30]);
      else if( id[0]==33&&id[1]==33)
	norm(_couplast*_coup[31]);
      else if((id[0]==32&&id[1]==33) || (id[0]==33&&id[1]==32))
	norm(_couplast*_coup[32]);
      else if((id[0]==22&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[33]);
      else if((id[0]==23&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[34]);
      else
	assert(false);
    }
  }
}
