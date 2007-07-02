// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMFFHVertex class.
//

#include "NMSSMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Models/Susy/NMSSM.h"

using namespace Herwig;
using namespace Herwig::Helicity;

NMSSMFFHVertex::NMSSMFFHVertex() {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks and neutral higgs
  int in[5]={25,35,45,36,46};
  for(unsigned int iy=0;iy<5;++iy) {
    for(unsigned int ix=1;ix<7;++ix) {
      first.push_back(-ix);
      second.push_back(ix);
      third.push_back(in[iy]);
    }
  }
  // leptons and neutral higgs
  for(unsigned int iy=0;iy<5;++iy) {
    for(unsigned int ix=11;ix<17;ix+=2) {
      first.push_back(-ix);
      second.push_back(ix);
      third.push_back(in[iy]);
    }
  }
  // the quarks  and the charged higgs
  for(unsigned int ix=0;ix<3;++ix) {
    first.push_back(2*ix);
    second.push_back(-2*ix-1);
    third.push_back(-37);
  }
  // the leptons and the charged higgs
  for(unsigned int ix=0;ix<3;++ix) {
    first.push_back(2*ix+12);
    second.push_back(-2*ix-11);
    third.push_back(-37);
  }
  setList(first,second,third);
  _mw=0.;
  _sinb=0.;
  _cosb=0.;
  _tanb=0.;
  _idlast=make_pair(0,0);
  _q2last=0.;
  _masslast=make_pair(0.,0);
  _couplast=0.;
}

void NMSSMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _mixS << _mixP << _mw << _sinb << _cosb << _tanb << _sw << _theSM;
}

void NMSSMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixS >> _mixP >> _mw >> _sinb >> _cosb >> _tanb >> _sw >> _theSM;
}

void NMSSMFFHVertex::doinit() throw(InitException) {
  FFSVertex::doinit();
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  _theSM = model;
  // sin theta_W
  double sw2=_theSM->sin2ThetaW();
  _sw = sqrt(sw2);
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  // Mass of the W boson
  _mw=getParticleData(ParticleID::Wplus)->mass();
  // sin and cos beta
  _tanb = model->tanBeta();
  double beta = atan(_tanb);
  _sinb=sin(beta);
  _cosb=cos(beta);
  cerr << *this << "\n";
}

ClassDescription<NMSSMFFHVertex> NMSSMFFHVertex::initNMSSMFFHVertex;
// Definition of the static class description member.

void NMSSMFFHVertex::Init() {

  static ClassDocumentation<NMSSMFFHVertex> documentation
    ("The NMSSMFFHVertex class implements the vertex for the couplings"
     " of the Higgs bosons of the NMSSM to Standard Model fermions");

}

void NMSSMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c, int) {
  int ihiggs=c->id();
  int id(abs(a->id()));
  Complex output(1.);
  // neutral Higgs
  if(ihiggs==25||ihiggs==35||ihiggs==45||ihiggs==36||ihiggs==46) {
    // left and right couplings set to one
    if(_idlast.first!=id||q2!=_q2last) {
      _idlast.first=id;
      _masslast.first = _theSM->mass(q2,a);
    }
    output*=_masslast.first;
    // CP-even 
    if(ihiggs==25||ihiggs==35||ihiggs==45) {
      int iloc = (ihiggs-25)/10;
      output *= (id%2==0) ? (*_mixS)(iloc,1)/_sinb : (*_mixS)(iloc,0)/_cosb;
      setLeft(1.); setRight(1.);
    } 
    // CP-odd
    else {
      int iloc = (ihiggs-36)/10;
      output *= (id%2==0) ? (*_mixP)(iloc,1)/_sinb : (*_mixP)(iloc,0)/_cosb;
      setLeft(1.); setRight(-1.);
      output *= Complex(0.,1.);
      Complex fact = (id%2==0) ? (*_mixP)(iloc,1)/_sinb : (*_mixP)(iloc,0)/_cosb;
    }
  }
  // charged higgs
  else if(abs(ihiggs)==37) {
    output*=-sqrt(2.);
    int id2=abs(b->id());
    if(id2<id) swap(id,id2);
    if(_idlast.first!=id||_idlast.second!=id2||q2!=_q2last) {
      _idlast.first =id ;
      _idlast.second=id2;
      _masslast.first  = _theSM->mass(q2,a);
      _masslast.second = _theSM->mass(q2,b);
    }
    if(ihiggs>0) {
      setRight(_masslast.first *_tanb);
      setLeft (_masslast.second/_tanb);
    }
    else {
      setLeft (_masslast.first *_tanb);
      setRight(_masslast.second/_tanb);
    }
  }
  else {
    throw Exception() << "Unknown Higgs boson, PDG code = " << ihiggs 
		      << "in NMSSMFFHVertex::setCoupling()"
		      << Exception::runerror;
  }
  // prefact
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = 0.5*sqrt(4.0*pi*alpha)/_sw/_mw;
    _q2last=q2;
  }
  setNorm(-_couplast*output);
}
