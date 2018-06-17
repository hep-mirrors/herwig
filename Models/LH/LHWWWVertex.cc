// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWWVertex class.
//

#include "LHWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _corr; 
}

void LHWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _corr; 
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWWVertex,VVVVertex>
describeHerwigLHWWWVertex("Herwig::LHWWWVertex", "HwLHModel.so");

void LHWWWVertex::Init() {

  static ClassDocumentation<LHWWWVertex> documentation
    ("The LHWWWVertex class implements the triple electroweak"
     " gauge boson couplings in the Little Higgs model.");

}

LHWWWVertex::LHWWWVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHWWWVertex::doinit() {
  // particles
  addToList(24,  -24,  22);
  addToList(24,  -24,  23);
  addToList(24,  -24,  32);
  addToList(24,  -24,  33);

  addToList(34,  -24,  23);
  addToList(34,  -24,  32);
  addToList(34,  -24,  33);

  addToList(24,  -34,  23);
  addToList(24,  -34,  32);
  addToList(24,  -34,  33);

  addToList(34,  -34,  22);
  addToList(34,  -34,  23);
  addToList(34,  -34,  32);
  addToList(34,  -34,  33);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWVertex::doinit()"
			  << Exception::runerror;
  // correction factors for the different interactions
  double sw(sqrt(model->sin2ThetaW())),cw(sqrt(1.-model->sin2ThetaW()));
  double vf(sqr(model->vev()/model->f()));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double xB(-2.5/sw*sp*cp*(sqr(cp)-sqr(sp)));
  double xH(2.5/sw/cw*s*c*sp*cp*(sqr(c*sp)+sqr(s*cp))/
	    (5.*sqr(sp*cp/sw)-sqr(s*c/cw)));
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  _corr.resize(12);
  // W_L W_L A_L
  _corr[ 0] = -1.;
  // W_L W_L A_H
  _corr[ 1] = cw/sw*vf*xB;
  // W_L W_H A_L
  _corr[ 2] = 0.;
  // W_L W_H A_H
  _corr[ 3] = -vf/sw*xH;
  // W_H W_H A_L
  _corr[ 4] = -1.;
  // W_H W_H A_H
  _corr[ 5] = vf/sw*(xH*(sqr(c)-sqr(s))/s/c+cw*xB);
  // W_L W_L Z_L
  _corr[ 6] = -cw/sw;
  // W_L W_L Z_H
  _corr[ 7] = vf/sw*(cw*xW+s*c*(sqr(c)-sqr(s)));
  // W_L W_H Z_L
  _corr[ 8] = -vf/sw*xW;
  // W_L W_H Z_H
  _corr[ 9] = -1./sw;
  // W_H W_H Z_L
  _corr[10] = -cw/sw;
  // W_H W_H Z_H
  _corr[11] = (sqr(c)-sqr(s))/s/c/sw;
  VVVVertex::doinit();
}

// couplings for the WWW vertex
void LHWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  int ia(a->iCharge()/3),ib(b->iCharge()/3),ic(c->iCharge()/3);
  int ida(a->id()), idb(b->id()), idc(c->id());
  // get the particles in the interaction
  int ineut,nh(0);
  if(ia==0) {
    ineut=ida;
    if(abs(idb)==34) ++nh;
    if(abs(idc)==34) ++nh;
  }
  else if(ib==0) {
    ineut=idb;
    if(abs(ida)==34) ++nh;
    if(abs(idc)==34) ++nh;
  }
  else {
    ineut=idc;
    if(abs(ida)==34) ++nh;
    if(abs(idb)==34) ++nh;
  }
  if(nh==0) {
    if     (ineut==22) norm(_corr[ 0]*_couplast);
    else if(ineut==23) norm(_corr[ 6]*_couplast);
    else if(ineut==32) norm(_corr[ 1]*_couplast);
    else if(ineut==33) norm(_corr[ 7]*_couplast);
  }
  else if(nh==1) {
    if     (ineut==22) norm(_corr[ 2]*_couplast);
    else if(ineut==23) norm(_corr[ 8]*_couplast);
    else if(ineut==32) norm(_corr[ 3]*_couplast);
    else if(ineut==33) norm(_corr[ 9]*_couplast);
  }
  else if(nh==2) {
    if     (ineut==22) norm(_corr[ 4]*_couplast);
    else if(ineut==23) norm(_corr[10]*_couplast);
    else if(ineut==32) norm(_corr[ 5]*_couplast);
    else if(ineut==33) norm(_corr[11]*_couplast);
  }
  // check the order for the overall sign 
  if((ia<0  && ib>0  && ic==0) || 
     (ia==0 && ib<0  && ic>0 ) || 
     (ia>0  && ib==0 && ic<0 ) ) norm(-norm());
}
