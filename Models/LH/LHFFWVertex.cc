// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFWVertex class.
//

#include "LHFFWVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

LHFFWVertex::LHFFWVertex() 
  : _ckm(3,vector<Complex>(3,0.0)), _couplast(0.), _q2last(0.*GeV2),
    _corrL(0.),_corrH(0.),_tcorrL(0.),_tcorrH(0.),_tHcorrL(0.), _tHcorrH(0.) {
  // order of vertex in couplings
  orderInGem(1);
  orderInGs(0);
}

void LHFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _ckm << _corrL << _corrH << _tcorrL << _tcorrH << _tHcorrL << _tHcorrH;
}

void LHFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _ckm >> _corrL >> _corrH >> _tcorrL >> _tcorrH >> _tHcorrL >> _tHcorrH;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFWVertex,FFVVertex>
describeHerwigLHFFWVertex("Herwig::LHFFWVertex", "HwLHModel.so");

void LHFFWVertex::Init() {

  static ClassDocumentation<LHFFWVertex> documentation
    ("The LHFFWVertex class implements the vertices for"
     " the coupling of the W and heavy W to the Standard Model "
     "fermions and the heavy top quark in the Little Higgs model");

}
  
void LHFFWVertex::doinit() {
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      addToList(-ix,      iy,      -24);
      addToList(-ix,      iy,      -34);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix,    ix+1,    -24);
    addToList(-ix,    ix+1,    -34);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      addToList(-ix,      iy,      24);
      addToList(-ix,      iy,      34);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1,    ix,    24);
    addToList(-ix-1,    ix,    34);
  }
  // couplings to new heavy quark
  addToList(-5,  8,  -24);
  addToList(-5,  8,  -34);
  addToList(-8,  5,  24);
  addToList(-8,  5,  34);
  ThePEG::Helicity::FFVVertex::doinit();
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHFFWVertex::doinit()"
			  << Exception::runerror;
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(generator()->standardModel()->CKM());
  if(!hwCKM) 
    throw InitException() << "Must have access to the Herwig::StandardCKM object"
			  << "for the CKM matrix in LHFFWVertex::doinit()"
			  << Exception::runerror;
  _ckm = hwCKM->getUnsquaredMatrix(model->families());
  // compute the correction factors
  double s2(sqr(model->sinTheta())),c2(sqr(model->cosTheta()));
  double vf(model->vev()/model->f());
  double xL = sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2()));
  // from Table VIII with -sign to agree with our SM conventions
  _corrL   =  1.-0.5*sqr(vf)*c2*(c2-s2);
  _corrH   = -model->cosTheta()/model->sinTheta();
  _tcorrL  =  1.-0.5*sqr(vf)*(c2*(c2-s2)+sqr(xL));
  _tcorrH  = -model->cosTheta()/model->sinTheta();
  _tHcorrL = -vf*xL;
  _tHcorrH =  vf*xL*model->cosTheta()/model->sinTheta();
}

void LHFFWVertex::setCoupling(Energy2 q2, tcPDPtr a, 
				       tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast    = -sqrt(0.5)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  right(0.);
  // the left and right couplings
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  bool heavy(false);
  if(iferm==8) {
    iferm = 6;
    heavy = true;
  }
  if(ianti==8) {
    ianti = 6;
    heavy = true;
  }
  // quarks
  if(iferm>=1 && iferm <=6) {
    int iu,id;
    // up type first
    if(iferm%2==0) {
      iu = iferm/2;
      id = (ianti+1)/2;
    }
    // down type first
    else {
      iu = ianti/2;
      id = (iferm+1)/2;
    }
    assert( iu>=1 && iu<=3 && id>=1 && id<=3);
    left(_ckm[iu-1][id-1]);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    left(1.);
  }
  else assert(false);
  // correction factors
  // light W
  if(abs(c->id())==ParticleID::Wplus) {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      left(_corrL*left());
    }
    // light top quark
    else if(!heavy) {
      left(_tcorrL*left());
    }
    // heavy top quark
    else {
      left(_tHcorrL*left());
    }
  }
  // heavy W
  else {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      left(_corrH*left());
    }
    // light top quark
    else if(!heavy) {
      left(_tcorrH*left());
    }
    // heavy top quark
    else {
      left(_tHcorrH*left());
    }
  }
}
