// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LittleHiggsFFHVertex class.
//

#include "LittleHiggsFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LittleHiggsFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _model << ounit(_coup,1./GeV);
}

void LittleHiggsFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _model >> iunit(_coup,1./GeV);
}

ClassDescription<LittleHiggsFFHVertex> LittleHiggsFFHVertex::initLittleHiggsFFHVertex;
// Definition of the static class description member.

void LittleHiggsFFHVertex::Init() {

  static ClassDocumentation<LittleHiggsFFHVertex> documentation
    ("The LittleHiggsFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model");

}

LittleHiggsFFHVertex::LittleHiggsFFHVertex() 
  : _idlast(0), _q2last(0.*GeV2), _masslast(0.*GeV) {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(35);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(35);
  }
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(36);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(36);
  }




  setList(first,second,third);
//   // set up for the couplings
//   _couplast=InvEnergy();
//   _idlast=0;
//   _q2last=0.*GeV2;
//   _masslast=0.*GeV;
//   _mw=0.*GeV;
//   _sw=0.;
}

void LittleHiggsFFHVertex::doinit() throw(InitException) {
//   _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
//   if (!_theSM) 
//     throw InitException();
//   double sw2=_theSM->sin2ThetaW();
//   _sw = sqrt(sw2);
//   _mw= getParticleData(ThePEG::ParticleID::Wplus)->mass();
//   orderInGem(1);
//   orderInGs(0);
//   FFSVertex::doinit();
}

void LittleHiggsFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr, int) {
//   int iferm=abs(a->id());
//   // left and right couplings set to one
//   setLeft(1.); setRight(1.);
//   // first the overall normalisation
//   if(q2!=_q2last) {
//     double alpha = _theSM->alphaEM(q2);
//     _couplast = -0.5*sqrt(4.0*Constants::pi*alpha)/_sw/_mw;
//     _q2last=q2;
//     _idlast=iferm;
//     if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
//       _masslast=_theSM->mass(q2,a);
//     }
//     else {
//       throw HelicityConsistencyError() << "LittleHiggsFFHVertex::setCoupling " 
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
//       throw HelicityConsistencyError() << "LittleHiggsFFHVertex::setCoupling " 
// 				       << "Unknown particle in Higgs vertex" 
// 				       << Exception::warning;
//       _masslast = 0*MeV;
//     }
//   }
//   setNorm(_couplast*_masslast);
}
