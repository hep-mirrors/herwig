// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
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

IBPtr FFSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFSDecayer::doinit() {
  _perturbativeVertex        = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex            = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  _abstractIncomingVertex    = dynamic_ptr_cast<AbstractFFVVertexPtr>(getIncomingVertex());

  if (getOutgoingVertices()[0]){
    if (getOutgoingVertices()[0]->getName()==VertexType::FFV){
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
      _abstractOutgoingVertexS   = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[1]);
    }
    else {
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
      _abstractOutgoingVertexS   = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[0]);
    }
  }
  else if (getOutgoingVertices()[1]){
    if (getOutgoingVertices()[1]->getName()==VertexType::FFV){
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[1]);
      _abstractOutgoingVertexS   = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[0]);
    }
    else {
      _abstractOutgoingVertexF   = dynamic_ptr_cast<AbstractFFVVertexPtr>(getOutgoingVertices()[0]);
      _abstractOutgoingVertexS   = dynamic_ptr_cast<AbstractVSSVertexPtr>(getOutgoingVertices()[1]);
    }
  }
  GeneralTwoBodyDecayer::doinit();
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _perturbativeVertex       << _abstractVertex
     << _abstractIncomingVertex   << _abstractOutgoingVertexF
     << _abstractOutgoingVertexS;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _perturbativeVertex       >> _abstractVertex
     >> _abstractIncomingVertex   >> _abstractOutgoingVertexF
     >> _abstractOutgoingVertexS;
}

ClassDescription<FFSDecayer> FFSDecayer::initFFSDecayer;
// Definition of the static class description member.

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

double FFSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id()    > 0? 0:1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC())     itype[1] = decay[0]->id() > 0? 0:1;
  else                              itype[1] = 2;
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));

  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(_wave,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wave[0].wave().Type() != u_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) _wave   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wavebar[0].wave().Type() != v_spinortype)
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
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
  }
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(_wavebar,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(_wave   ,decay[0],outgoing);
  ScalarWaveFunction scal(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) (*ME())(if1, if2, 0) = 
	_abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],scal);
      else     (*ME())(if2, if1, 0) = 
	_abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],scal);
    }
  }
  double output = (ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    double mu1(0.),mu2(0.);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, outa.first, outb.first);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, outb.first, outa.first);
      
    }
    double c2 = norm(_perturbativeVertex->norm());
    Complex cl = _perturbativeVertex->left();
    Complex cr = _perturbativeVertex->right();
    double me2 = c2*( (norm(cl) + norm(cr))*(1. + sqr(mu1) - sqr(mu2))
		      + 2.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/16./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double FFSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {

  int iscal (0), iferm (1), iglu (2);
  // get location of outgoing fermion/scalar
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,iferm);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(itype[0] == 0 || itype[1] == 0 || 
	   (itype[0] == 2 && itype[1] == 2 && decay[iscal]->id() < 0));  

  if(meopt==Initialize) {
    // create spinor (bar) for decaying particle
    if(ferm) {
      SpinorWaveFunction::calculateWaveFunctions(_wave3, _rho3, const_ptr_cast<tPPtr>(&inpart), 
						 incoming);
      if(_wave3[0].wave().Type() != u_spinortype)
   	for(unsigned int ix = 0; ix < 2; ++ix) _wave3[ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar3,_rho3, const_ptr_cast<tPPtr>(&inpart), 
						    incoming);
      if(_wavebar3[0].wave().Type() != v_spinortype)
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
    ScalarWaveFunction::constructSpinInfo(        decay[iscal],outgoing,true);
    VectorWaveFunction::constructSpinInfo(_gluon, decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1] = cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin0,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(_wavebar3, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(_wave3   , decay[iferm],outgoing);
  
  ScalarWaveFunction _swave3(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(_gluon,   decay[iglu ],outgoing,true);

  // // gauge invariance test
  //   _gluon.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) _gluon.push_back(VectorWaveFunction());
  //   else {
  //     _gluon.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 					  decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  if (! ((_abstractIncomingVertex  && (_abstractOutgoingVertexF || _abstractOutgoingVertexS)) ||
	 (_abstractOutgoingVertexF &&  _abstractOutgoingVertexS)))
    throw Exception()
      << "Invalid vertices for QCD radiation in FFS decay in FFSDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), S(2);
  if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 && 
      decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);
  else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
	   decay[iscal]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);


  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
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
		<< spinorInter.particle()->PDGName()  << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = _abstractVertex->evaluate(scale,spinorInter,_wavebar3[ifo],_swave3)/gs;
	  }
	  else {
	    SpinorBarWaveFunction spinorBarInter = 
	      _abstractIncomingVertex->evaluate(scale,3,inpart.dataPtr(),_wavebar3[ifi],
					       _gluon[2*ig],inpart.mass());

	    if (_wavebar3[ifi].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
	      throw Exception()
		<< _wavebar3[ifi].particle()->PDGName()  << " was changed to " 
		<< spinorBarInter.particle()->PDGName()  << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = _abstractVertex->evaluate(scale,_wave3[ifo], spinorBarInter,_swave3)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(ifi, 0, ifo, ig) += 
	       colourFlow[0][ix].second*diag;
	  }
	}
	  
  	// radiation from outgoing fermion
  	if(decay[iferm]->dataPtr()->coloured()) {
  	  assert(_abstractOutgoingVertexF);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();

	  double gs   = _abstractOutgoingVertexF->strongCoupling(scale);	  	  
	  if (ferm) {	    
	    SpinorBarWaveFunction spinorBarInter = 
	      _abstractOutgoingVertexF->evaluate(scale,3,off,_wavebar3[ifo],
						_gluon[2*ig],decay[iferm]->mass());
	    
	    if(_wavebar3[ifo].particle()->PDGName()!=spinorBarInter.particle()->PDGName())
	      throw Exception()
		<< _wavebar3[ifo].particle()->PDGName() << " was changed to " 
		<< spinorBarInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = _abstractVertex->evaluate(scale,_wave3[ifi],spinorBarInter,_swave3)/gs;
	  }
	  else {
	    SpinorWaveFunction spinorInter = 
	      _abstractOutgoingVertexF->evaluate(scale,3,off,_wave3[ifo],
						_gluon[2*ig],decay[iferm]->mass());
	      
	    if(_wave3[ifo].particle()->PDGName()!=spinorInter.particle()->PDGName())
	      throw Exception()
		<< _wave3[ifo].particle()->PDGName() << " was changed to " 
		<< spinorInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
		<< Exception::runerror;
	    diag = _abstractVertex->evaluate(scale,spinorInter,_wavebar3[ifi],_swave3)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	    (*ME[colourFlow[F][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[F][ix].second*diag;
	  }
  	}

  	// radiation from outgoing scalar
  	if(decay[iscal]->dataPtr()->coloured()) {
  	  assert(_abstractOutgoingVertexS);
	  // ensure you get correct ougoing particle from first vertex
	  tcPDPtr off = decay[iscal]->dataPtr();
	  if(off->CC()) off = off->CC();
	  
	  double gs = _abstractOutgoingVertexS->strongCoupling(scale);
	  ScalarWaveFunction  scalarInter = 
	    _abstractOutgoingVertexS->evaluate(scale,3,off,_gluon[2*ig],
					      _swave3,decay[iscal]->mass());
	    
	  if(_swave3.particle()->PDGName()!=scalarInter.particle()->PDGName())
	    throw Exception()
	      << _swave3    .particle()->PDGName() << " was changed to " 
	      << scalarInter.particle()->PDGName() << " in FFSDecayer::threeBodyME"
	      << Exception::runerror; 
	  if (ferm){
	    diag = _abstractVertex->evaluate(scale,_wave3[ifi],_wavebar3[ifo],scalarInter)/gs;
	  }
	  else {
	    diag = _abstractVertex->evaluate(scale,_wave3[ifo],_wavebar3[ifi],scalarInter)/gs;
	  }
	  for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
  	    (*ME[colourFlow[S][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[S][ix].second*diag;
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
