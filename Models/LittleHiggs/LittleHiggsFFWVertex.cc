// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsFFWVertex class.
//

#include "LittleHiggsFFWVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

inline LittleHiggsFFWVertex::LittleHiggsFFWVertex() 
  : _ckm(3,vector<Complex>(3,0.0)), _couplast(0.), _q2last(0.*GeV2),
    _corrL(0.),_corrH(0.),_tcorrL(0.),_tcorrH(0.),_tHcorrL(0.), _tHcorrH(0.) {
  // particles for the vertex
  vector<int> first,second,third;
  // particles for outgoing W-
  // quarks
  for(unsigned int ix=1;ix<6;ix+=2) {
    for(unsigned int iy=2;iy<7;iy+=2) {
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(-24);
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(-34);
    }
  }
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix+1);
    third.push_back(-24);
    first.push_back(-ix);
    second.push_back(ix+1);
    third.push_back(-34);
  }
  // particles for outgoing W+
  // quarks
  for(unsigned int ix=2;ix<7;ix+=2) {
    for(unsigned int iy=1;iy<6;iy+=2) {
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(24);
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(34);
    }
  }
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix-1);
    second.push_back(ix);
    third.push_back(24);
    first.push_back(-ix-1);
    second.push_back(ix);
    third.push_back(34);
  }
  // couplings to new heavy quark
  first.push_back(-5);
  second.push_back(8);
  third.push_back(-24);
  first.push_back(-5);
  second.push_back(8);
  third.push_back(-34);
  first.push_back(-8);
  second.push_back(5);
  third.push_back(24);
  first.push_back(-8);
  second.push_back(5);
  third.push_back(34);
  setList(first,second,third);
}

void LittleHiggsFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _model << _ckm << _corrL << _corrH << _tcorrL << _tcorrH 
     << _tHcorrL << _tHcorrH;
}

void LittleHiggsFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _ckm >> _corrL >> _corrH >> _tcorrL >> _tcorrH 
     >> _tHcorrL >> _tHcorrH;
}

ClassDescription<LittleHiggsFFWVertex> LittleHiggsFFWVertex::initLittleHiggsFFWVertex;
// Definition of the static class description member.

void LittleHiggsFFWVertex::Init() {

  static ClassDocumentation<LittleHiggsFFWVertex> documentation
    ("The LittleHiggsFFWVertex class implements the vertices for"
     " the coupling of the W and heavy W to the Standard Model "
     "fermions and the heavy top quark in the Little Higgs model");

}
  
void LittleHiggsFFWVertex::doinit() throw(InitException) {
  ThePEG::Helicity::FFVVertex::doinit();
  _model = dynamic_ptr_cast<cLittleHiggsModelPtr>(generator()->standardModel());
  if(!_model) 
    throw InitException() << "Must be using the LittleHiggsModel "
			  << " in LittleHiggsFFWVertex::doinit()"
			  << Exception::runerror;
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(generator()->standardModel()->CKM());
  if(!hwCKM) 
    throw InitException() << "Must have access to the Herwig::StandardCKM object"
			  << "for the CKM matrix in LittleHiggsFFWVertex::doinit()"
			  << Exception::runerror;
  _ckm = hwCKM->getUnsquaredMatrix(_model->families());
  // compute the correction factors
  double s2(sqr(_model->sinTheta())),c2(sqr(_model->cosTheta()));
  double vf(_model->vev()/_model->f());
  double xL = sqr(_model->lambda1())/(sqr(_model->lambda1())+sqr(_model->lambda2()));
  _corrL   =  1.-0.5*sqr(vf)*c2*(c2-s2);
  _corrH   = -_model->cosTheta()/_model->sinTheta();
  _tcorrL  =  1.-0.5*sqr(vf)*(c2*(c2-s2)+sqr(xL));
  _tcorrH  = -_model->cosTheta()/_model->sinTheta();
  _tHcorrL =  vf*xL;
  _tHcorrH = -vf*xL*_model->cosTheta()/_model->sinTheta();
  // order of vertex in couplings
  orderInGem(1);
  orderInGs(0);
}

void LittleHiggsFFWVertex::setCoupling(Energy2 q2, tcPDPtr a, 
				       tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _model->alphaEM(q2);
    double sw    = sqrt(2.*(_model->sin2ThetaW()));
    _couplast    = sqrt(4.0*Constants::pi*alpha)/sw;
    _q2last=q2;
  }
  setNorm(_couplast);
  setRight(0.);
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
    if( iu<1 || iu>3 || id<1 || id>3)
      throw HelicityConsistencyError() << "LittleHiggsFFWVertex::setCoupling "
				       << "Unknown particle in W vertex" 
				       << Exception::runerror;
    setLeft(_ckm[iu-1][id-1]);
    setRight(0.);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    setLeft(1.);
  }
  else
    throw HelicityConsistencyError() << "LittleHiggsFFWVertex::setCoupling "
				     << "Unknown particle in W vertex" 
				     << Exception::runerror;
  // correction factors
  // light W
  if(abs(c->id())==ParticleID::Wplus) {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      setLeft(_corrL*getLeft());
    }
    // light top quark
    else if(!heavy) {
      setLeft(_tcorrL*getLeft());
    }
    // heavy top quark
    else {
      setLeft(_tHcorrL*getLeft());
    }
  }
  // heavy W
  else {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      setLeft(_corrH*getLeft());
    }
    // light top quark
    else if(!heavy) {
      setLeft(_tcorrH*getLeft());
    }
    // heavy top quark
    else {
      setLeft(_tHcorrH*getLeft());
    }
  }
}
