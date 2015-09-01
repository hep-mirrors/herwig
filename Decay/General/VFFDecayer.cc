// -*- C++ -*-
//
// VFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VFFDecayer class.
//

#include "VFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void VFFDecayer::doinit() {
  _perturbativeVertex      = dynamic_ptr_cast<FFVVertexPtr>        (getVertex());
  _abstractVertex          = dynamic_ptr_cast<AbstractFFVVertexPtr>(getVertex());
  _abstractIncomingVertex  = dynamic_ptr_cast<AbstractVVVVertexPtr>(getIncomingVertex());
  _abstractOutgoingVertex1 = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
  _abstractOutgoingVertex2 = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
  GeneralTwoBodyDecayer::doinit();
}

void VFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex           << _perturbativeVertex
     << _abstractIncomingVertex   << _abstractOutgoingVertex1
     << _abstractOutgoingVertex2;
}

void VFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex           >> _perturbativeVertex
     >> _abstractIncomingVertex   >> _abstractOutgoingVertex1
     >> _abstractOutgoingVertex2;
}

ClassDescription<VFFDecayer> VFFDecayer::initVFFDecayer;
// Definition of the static class description member.

void VFFDecayer::Init() {

  static ClassDocumentation<VFFDecayer> documentation
    ("The VFFDecayer implements the matrix element for the"
     " decay of a vector to fermion-antifermion pair");

}

double VFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
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
  // compute the matrix element
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int ifm = 0; ifm < 2; ++ifm) { //loop over fermion helicities
    for(unsigned int ia = 0; ia < 2; ++ia) {// loop over antifermion helicities
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {//loop over vector helicities
	if(iferm > ianti) {
	  (*ME())(vhel, ia, ifm) = 
	    _abstractVertex->evaluate(scale,_wave[ia],
				      _wavebar[ifm],_vectors[vhel]);
	}
	else
	  (*ME())(vhel,ifm,ia)=
	    _abstractVertex->evaluate(scale,_wave[ia],
				      _wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    double mu1(outa.second/inpart.second), mu2(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(sqr(inpart.second), outa.first, outb.first,in);
    Complex cl(_perturbativeVertex->left()), cr(_perturbativeVertex->right());
    double me2 = (norm(cl) + norm(cr))*( sqr(sqr(mu1) - sqr(mu2)) 
					 + sqr(mu1) + sqr(mu2) - 2.)
      - 6.*(cl*conj(cr) + cr*conj(cl)).real()*mu1*mu2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(_perturbativeVertex->norm())*me2*pcm / 
      (24.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double VFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {

  bool massless = inpart.mass()==ZERO;

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
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(_vector3, _rho3, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, massless);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(_vector3 ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,massless);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar3,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave3   ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(_gluon   ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1]=cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,     PDT::Spin1Half,
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

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ifm = 0; ifm < 2; ++ifm) {
      for(unsigned int ia = 0; ia < 2; ++ia) {
	for(unsigned int ig = 0; ig < 2; ++ig) {
	  // radiation from the incoming vector
	  if(inpart.dataPtr()->coloured()) {
	    assert(_abstractIncomingVertex);
	    
	    VectorWaveFunction vectorInter = 
	      _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),_vector3[iv],
						_gluon[2*ig],inpart.mass());
	    
	    if (_vector3[iv].particle()->PDGName()!=vectorInter.particle()->PDGName())
	      throw Exception()
		<< _vector3[iv].particle()->PDGName() << " was changed to " 
		<< vectorInter .particle()->PDGName() << " in VFFDecayer::threeBodyME"
		<< Exception::runerror;

	    double gs    = _abstractIncomingVertex->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,_wave3[ia],_wavebar3[ifm],vectorInter)/gs;
	    for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	      (*ME[colourFlow[0][ix].first])(iv, ia, ifm, ig) += 
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
		<< interS        .particle()->PDGName() << " in VFFDecayer::threeBodyME"
		<< Exception::runerror;
	    
	    double gs    =  abstractOutgoingVertexF->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,_wave3[ia], interS,_vector3[iv])/gs;
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(iv, ia, ifm, ig) += 
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
		<< interS    .particle()->PDGName() << " in VFFDecayer::threeBodyME"
		<< Exception::runerror;

	    double gs    =  abstractOutgoingVertexA->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,interS,_wavebar3[ifm],_vector3[iv])/gs;
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(iv, ia, ifm, ig) += 
		colourFlow[2][ix].second*diag;
	    }
	  }
	}
      }
    }
    if(massless) ++iv;
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],_rho3)).real();
    }
  }
  output*=(4.*Constants::pi);

  //return output
  return output;
}


void VFFDecayer::identifyVertices(const int iferm, const int ianti,
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
    << "Invalid vertices for QCD radiation in VFF decay in VFFDecayer::identifyVertices"
    << Exception::runerror;
  
}

