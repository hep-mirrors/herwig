// -*- C++ -*-
//
// FFVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVDecayer class.
//

#include "FFVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVDecayer::doinit() {
  _perturbativeVertex     = dynamic_ptr_cast<FFVVertexPtr>        (getVertex());
  _abstractVertex         = dynamic_ptr_cast<AbstractFFVVertexPtr>(getVertex());
  _abstractIncomingVertex = dynamic_ptr_cast<AbstractFFVVertexPtr>(getIncomingVertex());

  if (getOutgoingVertices()[0]){
    if (getOutgoingVertices()[0]->getName()==VertexType::FFV){
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
      _abstractOutgoingVertexV   = dynamic_ptr_cast<AbstractVVVVertexPtr>(getOutgoingVertices()[1]);
    }
    else {
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
      _abstractOutgoingVertexV   = dynamic_ptr_cast<AbstractVVVVertexPtr>(getOutgoingVertices()[0]);
    }
  }
  else if (getOutgoingVertices()[1]){
    if (getOutgoingVertices()[1]->getName()==VertexType::FFV){
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
      _abstractOutgoingVertexV   = dynamic_ptr_cast<AbstractVVVVertexPtr>(getOutgoingVertices()[0]);
    }
    else {
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
      _abstractOutgoingVertexV   = dynamic_ptr_cast<AbstractVVVVertexPtr>(getOutgoingVertices()[1]);
    }
  }

  GeneralTwoBodyDecayer::doinit();
}

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex           << _perturbativeVertex
     << _abstractIncomingVertex   << _abstractOutgoingVertexF
     << _abstractOutgoingVertexV;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex           >> _perturbativeVertex
     >> _abstractIncomingVertex   >> _abstractOutgoingVertexF
     >> _abstractOutgoingVertexV;
}

double FFVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // type of process
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;  
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(_wave,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wave[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) _wave   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wavebar[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) _wavebar[ix].conjugate();
    }
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(_wave,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_wavebar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_wavebar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_wave,decay[0],outgoing,true);
    }
    VectorWaveFunction::
      constructSpinInfo(_vector,decay[1],outgoing,true,false);
  }
  Energy2 scale(sqr(inpart.mass()));
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(_wavebar,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(_wave   ,decay[0],outgoing);
  bool massless = decay[1]->dataPtr()->mass()==ZERO;
  VectorWaveFunction::
    calculateWaveFunctions(_vector,decay[1],outgoing,massless);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(massless && vhel == 1) ++vhel;
	if(ferm)
	  (*ME())(if1, if2,vhel) = 
	    _abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],_vector[vhel]);
	else
	  (*ME())(if2, if1, vhel) = 
	    _abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],_vector[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outa.first->iSpin() == PDT::Spin1Half)
      _perturbativeVertex->setCoupling(sqr(inpart.second), in,
				       outa.first, outb.first);
    else {
      swap(mu1,mu2);
      _perturbativeVertex->setCoupling(sqr(inpart.second),in,
				       outb.first,outa.first);
    }
    Complex cl(_perturbativeVertex->left()),cr(_perturbativeVertex->right());
    double me2(0.);
    if( mu2 > 0. ) {
      me2 = (norm(cl) + norm(cr))*(1. + sqr(mu1*mu2) + sqr(mu2) 
				   - 2.*sqr(mu1) - 2.*sqr(mu2*mu2) 
				   +  sqr(mu1*mu1))
	- 6.*mu1*sqr(mu2)*(conj(cl)*cr + conj(cr)*cl).real();
      me2 /= sqr(mu2);
    }
    else
      me2 = 2.*( (norm(cl) + norm(cr))*(sqr(mu1) + 1.) 
		 - 4.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = norm(_perturbativeVertex->norm())*me2*pcm/16./Constants::pi; 
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer 
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

ClassDescription<FFVDecayer> FFVDecayer::initFFVDecayer;
// Definition of the static class description member.

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("There is no documentation for the FFVDecayer class");

}

double  FFVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay, MEOption meopt) {

  int iferm (0), ivect (1), iglu (2);
  // get location of outgoing lepton/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin1Half) swap(iferm,ivect);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(itype[0] == 0 || itype[1] == 0 || 
	   (itype[0] == 2 && itype[1] == 2 && decay[ivect]->id() < 0));  

  // no emissions from massive vectors
  bool massless = decay[ivect]->dataPtr()->mass()==ZERO;
  if (_abstractOutgoingVertexV && (! massless))
    throw Exception()
      << "No dipoles available for massive vectors in FFVDecayer::threeBodyME"
      << Exception::runerror;

  if(meopt==Initialize) {
    // create spinor (bar) for decaying particle
    if(ferm) {
      SpinorWaveFunction::calculateWaveFunctions(_wave3, _rho3, const_ptr_cast<tPPtr>(&inpart), 
						 incoming);
      if(_wave3[0].wave().Type() != SpinorType::u)
   	for(unsigned int ix = 0; ix < 2; ++ix) _wave3[ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar3,_rho3, const_ptr_cast<tPPtr>(&inpart), 
						    incoming);
      if(_wavebar3[0].wave().Type() != SpinorType::v)
   	for(unsigned int ix = 0; ix < 2; ++ix) _wavebar3[ix].conjugate();
    }
  }
  // setup spin information when needed 
  if(meopt==Terminate) {
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(_wave3,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_wavebar3,decay[iferm],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_wavebar3,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_wave3,decay[iferm],outgoing,true);
    }
    VectorWaveFunction::constructSpinInfo(_vector3, decay[ivect],outgoing,true,massless);
    VectorWaveFunction::constructSpinInfo(_gluon,   decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1] = cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half,
								       PDT::Spin1,     PDT::Spin1)));

  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(_wavebar3, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(_wave3   , decay[iferm],outgoing);
  
  VectorWaveFunction::calculateWaveFunctions(_vector3, decay[ivect],outgoing,massless);
  VectorWaveFunction::calculateWaveFunctions(_gluon,   decay[iglu ],outgoing,true );

  // // gauge invariance test
  // _gluon.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) _gluon.push_back(VectorWaveFunction());
  //   else {
  //     _gluon.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 					  decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  if (! ((_abstractIncomingVertex  && (_abstractOutgoingVertexF  || _abstractOutgoingVertexV)) ||
	 (_abstractOutgoingVertexF &&  _abstractOutgoingVertexV)))
    throw Exception()
      << "Invalid vertices for QCD radiation in FFV decay in FFVDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), V(2);
  if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
      decay[ivect]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);
  else if (decay[ivect]->dataPtr()->iColour()==PDT::Colour3 && 
	   decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);

  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int iv = 0; iv < 3; ++iv) {
	for(unsigned int ig = 0; ig < 2; ++ig) {
	  // radiation from the incoming fermion
	  if(inpart.dataPtr()->coloured()) {
	    assert(_abstractIncomingVertex);
	    double gs = _abstractIncomingVertex->strongCoupling(scale);	  
	    if (ferm){
	      SpinorWaveFunction spinorInter =
		_abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),_wave3[ifi],
						  _gluon[2*ig],inpart.mass());
	      
	      if (_wave3[ifi].particle()->PDGName()!=spinorInter.particle()->PDGName())
		throw Exception()
		  << _wave3[ifi].particle()->PDGName()  << " was changed to " 
		  << spinorInter.particle()->PDGName()  << " in FFVDecayer::threeBodyME"
		  << Exception::runerror;
	      diag = _abstractVertex->evaluate(scale,spinorInter,_wavebar3[ifo],_vector3[iv])/gs;
	    }
	    else {
	      SpinorBarWaveFunction spinorBarInter = 
		_abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),_wavebar3[ifi],
						  _gluon[2*ig],inpart.mass());
	      
	      if (_wavebar3[ifi].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
		throw Exception()
		  << _wavebar3[ifi].particle()->PDGName()  << " was changed to " 
		  << spinorBarInter.particle()->PDGName()  << " in FFVDecayer::threeBodyME"
		  << Exception::runerror;
	      diag = _abstractVertex->evaluate(scale,_wave3[ifo], spinorBarInter,_vector3[iv])/gs;
	    }
	    for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	      (*ME[colourFlow[0][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[0][ix].second*diag;
	    }
	  }

	  // radiation from outgoing fermion
	  if(decay[iferm]->dataPtr()->coloured()) {
	    assert(_abstractOutgoingVertexF);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC();
	    
	    double gs = _abstractOutgoingVertexF->strongCoupling(scale);	  	  
	    if (ferm) {	    
	      SpinorBarWaveFunction spinorBarInter = 
		_abstractOutgoingVertexF->evaluate(scale,3,off,_wavebar3[ifo],
						   _gluon[2*ig],decay[iferm]->mass());
	      
	      if(_wavebar3[ifo].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
		throw Exception()
		  << _wavebar3[ifo].particle()->PDGName() << " was changed to " 
		  << spinorBarInter.particle()->PDGName() << " in FFVDecayer::threeBodyME"
		  << Exception::runerror;
	      diag = _abstractVertex->evaluate(scale,_wave3[ifi],spinorBarInter,_vector3[iv])/gs;
	    }
	    else {
	      SpinorWaveFunction spinorInter = 
		_abstractOutgoingVertexF->evaluate(scale,3,off,_wave3[ifo],
						   _gluon[2*ig],decay[iferm]->mass());
		
	      if(_wave3[ifo].particle()->PDGName()!=spinorInter.particle()->PDGName())
		throw Exception()
		  << _wave3[ifo].particle()->PDGName() << " was changed to " 
		  << spinorInter.particle()->PDGName() << " in FFVDecayer::threeBodyME"
		  << Exception::runerror;
	      diag = _abstractVertex->evaluate(scale,spinorInter,_wavebar3[ifi],_vector3[iv])/gs;
	    }
	    for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	      (*ME[colourFlow[F][ix].first])(ifi, ifo, iv, ig) += 
		 colourFlow[F][ix].second*diag;
	    }
	  }
	  
	  // radiation from outgoing vector
	  if(decay[ivect]->dataPtr()->coloured()) {
	    assert(_abstractOutgoingVertexV);
	    // ensure you get correct ougoing particle from first vertex
	    tcPDPtr off = decay[ivect]->dataPtr();
	    if(off->CC()) off = off->CC();

	    double sign  = decay[iferm]->id()>0 ? -1:1;
	    double gs    = _abstractOutgoingVertexV->strongCoupling(scale);
	    VectorWaveFunction  vectorInter = 
	      _abstractOutgoingVertexV->evaluate(scale,3,off,_gluon[2*ig],
						 _vector3[iv],decay[ivect]->mass());
	    
	    if(_vector3[iv].particle()->PDGName()!=vectorInter.particle()->PDGName())
	      throw Exception()
		<< _vector3[iv].particle()->PDGName() << " was changed to " 
		<< vectorInter. particle()->PDGName() << " in FFVDecayer::threeBodyME"
		<< Exception::runerror; 
	    if (ferm){
	      diag = sign*_abstractVertex->evaluate(scale,_wave3[ifi],_wavebar3[ifo],vectorInter)/gs;
	    }
	    else {
	      diag = sign*_abstractVertex->evaluate(scale,_wave3[ifo],_wavebar3[ifi],vectorInter)/gs;
	    }
	    for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	      (*ME[colourFlow[V][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[V][ix].second*diag;
	    }
	  }
	}
	if(massless) ++iv;
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
