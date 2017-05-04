// -*- C++ -*-
//
// SFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SFFDecayer::doinit() {
  _perturbativeVertex        = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex            = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  _abstractIncomingVertex    = dynamic_ptr_cast<AbstractVSSVertexPtr>(getIncomingVertex());
  _abstractOutgoingVertex1   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
  _abstractOutgoingVertex2   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
  GeneralTwoBodyDecayer::doinit();
}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex           << _perturbativeVertex 
     << _abstractIncomingVertex   << _abstractOutgoingVertex1
     << _abstractOutgoingVertex2;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex           >> _perturbativeVertex 
     >> _abstractIncomingVertex   >> _abstractOutgoingVertex1
     >> _abstractOutgoingVertex2;
}

ClassDescription<SFFDecayer> SFFDecayer::initSFFDecayer;
// Definition of the static class description member.

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

double SFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // work out which is the fermion and antifermion
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);

  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti){
	(*ME())(0, ia, ifm) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);
      }
      else {
	(*ME())(0, ifm, ia) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);	
      }
    }
  }

  double output = (ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(sqr(inpart.second), outb.first, outa.first,
				     in);
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    double c2 = norm(_perturbativeVertex->norm());
    Complex al(_perturbativeVertex->left()), ar(_perturbativeVertex->right());
    double me2 = -c2*( (norm(al) + norm(ar))*( sqr(mu1) + sqr(mu2) - 1.)
		       + 2.*(ar*conj(al) + al*conj(ar)).real()*mu1*mu2 );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/(8*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {
  
  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(iferm, ianti);
  if(itype[0]==2 && itype[1]==1) swap(iferm, ianti);
  if(itype[0]==0 && itype[1]==0 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  if(itype[0]==1 && itype[1]==1 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  
  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::
      calculateWaveFunctions(_rho3,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave3 = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar3 ,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave3    ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(_gluon    ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1]=cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,     PDT::Spin1Half,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar3, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave3   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(_gluon   , decay[iglu ],outgoing,true);

  // // gauge invariance test
  // _gluon.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) _gluon.push_back(VectorWaveFunction());
  //   else {
  //     _gluon.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 				          decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr abstractOutgoingVertexF;
  AbstractFFVVertexPtr abstractOutgoingVertexA;
  identifyVertices(iferm, ianti, inpart, decay, abstractOutgoingVertexF, abstractOutgoingVertexA);

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm) {
    for(unsigned int ia = 0; ia < 2; ++ia) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
	// radiation from the incoming scalar
	if(inpart.dataPtr()->coloured()) {
	  assert(_abstractIncomingVertex);

	  ScalarWaveFunction scalarInter = 
	    _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),
					      _gluon[2*ig],_swave3,inpart.mass());

	  if (_swave3.particle()->PDGName()!=scalarInter.particle()->PDGName())
	    throw Exception()
	      << _swave3    .particle()->PDGName() << " was changed to " 
	      << scalarInter.particle()->PDGName() << " in SFFDecayer::threeBodyME"
	      << Exception::runerror;

	  double gs    = _abstractIncomingVertex->strongCoupling(scale);
	  Complex diag = _abstractVertex->evaluate(scale,_wave3[ia], _wavebar3[ifm],
							  scalarInter)/gs;
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(0, ia, ifm, ig) += 
	       colourFlow[0][ix].second*diag;
	  }
	}	

	// radiation from outgoing fermion
	if(decay[iferm]->dataPtr()->coloured()) {
	  assert(abstractOutgoingVertexF);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorBarWaveFunction interS = 
	    abstractOutgoingVertexF->evaluate(scale,3,off,_wavebar3[ifm],
					       _gluon[2*ig],decay[iferm]->mass());
	  
	  if(_wavebar3[ifm].particle()->PDGName()!=interS.particle()->PDGName())
	    throw Exception()
	      << _wavebar3[ifm].particle()->PDGName() << " was changed to " 
	      << interS        .particle()->PDGName() << " in SFFDecayer::threeBodyME"
	      << Exception::runerror;

	  double gs    =  abstractOutgoingVertexF->strongCoupling(scale);
	  Complex diag = _abstractVertex->evaluate(scale,_wave3[ia], interS,_swave3)/gs;
	  for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	    (*ME[colourFlow[1][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[1][ix].second*diag;
	  }
	}

	// radiation from outgoing antifermion
	if(decay[ianti]->dataPtr()->coloured()) {
	  assert(abstractOutgoingVertexA);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[ianti]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorWaveFunction  interS = 
	    abstractOutgoingVertexA->evaluate(scale,3,off,_wave3[ia],
					      _gluon[2*ig],decay[ianti]->mass());

	  if(_wave3[ia].particle()->PDGName()!=interS.particle()->PDGName())
	    throw Exception()
	      << _wave3[ia].particle()->PDGName() << " was changed to " 
	      << interS    .particle()->PDGName() << " in SFFDecayer::threeBodyME"
	      << Exception::runerror;

	  double gs    =  abstractOutgoingVertexA->strongCoupling(scale);
	  Complex diag = _abstractVertex->evaluate(scale,interS,_wavebar3[ifm],_swave3)/gs;
	  for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	    (*ME[colourFlow[2][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[2][ix].second*diag;
	  }
	}
      }
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],_rho3)).real();
    }
  }
  output*=(4.*Constants::pi);

  // return the answer
  return output;
}

void SFFDecayer::identifyVertices(const int iferm, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractFFVVertexPtr & abstractOutgoingVertexF, 
				  AbstractFFVVertexPtr & abstractOutgoingVertexA){

  // work out which fermion each outgoing vertex corresponds to 
  // two outgoing vertices
  if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
    ((decay[iferm]->dataPtr()->iColour()==PDT::Colour3     &&
      decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
     (decay[iferm]->dataPtr()->iColour()==PDT::Colour8     &&
      decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
    if(_abstractOutgoingVertex1==_abstractOutgoingVertex2){
      abstractOutgoingVertexF = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iferm]->id()))){
      abstractOutgoingVertexF = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex2->isIncoming(getParticleData(decay[iferm]->id()))){
      abstractOutgoingVertexF = _abstractOutgoingVertex2;
      abstractOutgoingVertexA = _abstractOutgoingVertex1;
    }
  }
  else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	  decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
    if(_abstractOutgoingVertex1==_abstractOutgoingVertex2){
      abstractOutgoingVertexF = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iferm]->id()))){
      abstractOutgoingVertexF = _abstractOutgoingVertex1;
      abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (_abstractOutgoingVertex2->isIncoming(getParticleData(decay[iferm]->id()))){
      abstractOutgoingVertexF = _abstractOutgoingVertex2;
      abstractOutgoingVertexA = _abstractOutgoingVertex1;
    }
  }

  // one outgoing vertex
  else if(inpart.dataPtr()->iColour()==PDT::Colour3){
    if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
       decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
      if     (_abstractOutgoingVertex1) abstractOutgoingVertexF = _abstractOutgoingVertex1;
      else if(_abstractOutgoingVertex2) abstractOutgoingVertexF = _abstractOutgoingVertex2;
    }
    else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	     decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
      if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[ianti]->dataPtr()->id()))){
	abstractOutgoingVertexF = _abstractOutgoingVertex2;
	abstractOutgoingVertexA = _abstractOutgoingVertex1;
      }
      else {
	abstractOutgoingVertexF = _abstractOutgoingVertex1;
	abstractOutgoingVertexA = _abstractOutgoingVertex2;
      }
    }
  }
  else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
    if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
       decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
      if     (_abstractOutgoingVertex1) abstractOutgoingVertexA = _abstractOutgoingVertex1;
      else if(_abstractOutgoingVertex2) abstractOutgoingVertexA = _abstractOutgoingVertex2;
    }
    else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour8 &&
	     decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if (_abstractOutgoingVertex1->isIncoming(getParticleData(decay[iferm]->dataPtr()->id()))){
	abstractOutgoingVertexF = _abstractOutgoingVertex1;
	abstractOutgoingVertexA = _abstractOutgoingVertex2;
      }
      else {
	abstractOutgoingVertexF = _abstractOutgoingVertex2;
	abstractOutgoingVertexA = _abstractOutgoingVertex1;
      }
    }
  }
  
  if (! ((_abstractIncomingVertex  && (abstractOutgoingVertexF  || abstractOutgoingVertexA)) ||
	 ( abstractOutgoingVertexF &&  abstractOutgoingVertexA)))
    throw Exception()
    << "Invalid vertices for QCD radiation in SFF decay in SFFDecayer::identifyVertices"
    << Exception::runerror;

}
