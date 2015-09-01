// -*- C++ -*-
//
// TFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TFFDecayer class.
//

#include "TFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void TFFDecayer::doinit() {
  _perturbativeVertex        = dynamic_ptr_cast<FFTVertexPtr>         (getVertex());
  _abstractVertex            = dynamic_ptr_cast<AbstractFFTVertexPtr> (getVertex());
  _abstractOutgoingVertex1   = dynamic_ptr_cast<AbstractFFVVertexPtr> (getOutgoingVertices()[0]);
  _abstractOutgoingVertex2   = dynamic_ptr_cast<AbstractFFVVertexPtr> (getOutgoingVertices()[1]);
  _abstractFourPointVertex   = dynamic_ptr_cast<AbstractFFVTVertexPtr>(getFourPointVertex());

  GeneralTwoBodyDecayer::doinit();
}

void TFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex          << _perturbativeVertex
     << _abstractOutgoingVertex1 << _abstractOutgoingVertex2
     << _abstractFourPointVertex;
}

void TFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex          >> _perturbativeVertex
     >> _abstractOutgoingVertex1 >> _abstractOutgoingVertex2
     >> _abstractFourPointVertex;
}

ClassDescription<TFFDecayer> TFFDecayer::initTFFDecayer;
// Definition of the static class description member.

void TFFDecayer::Init() {

  static ClassDocumentation<TFFDecayer> documentation
    ("The TFFDecayer class implements the decay of a tensor particle "
     "to 2 fermions ");
  
}

double TFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  unsigned int iferm(0),ianti(1);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half)));
  if(decay[0]->id()>=0) swap(iferm,ianti);
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,_rho,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
  }
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(_tensors,const_ptr_cast<tPPtr>(&inpart),
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
  Energy2 scale(sqr(inpart.mass()));
  unsigned int thel,fhel,ahel;
  for(thel=0;thel<5;++thel) {
    for(fhel=0;fhel<2;++fhel) {
      for(ahel=0;ahel<2;++ahel) {
	if(iferm > ianti) {
	  (*ME())(thel,fhel,ahel) = 
	    _abstractVertex->evaluate(scale,_wave[ahel],
				      _wavebar[fhel],_tensors[thel]);
	}
	else {
	  (*ME())(thel,ahel,fhel) = 
	    _abstractVertex->evaluate(scale,_wave[ahel],
				      _wavebar[fhel],_tensors[thel]);
	}
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

Energy TFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale = sqr(inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(scale, in, outa.first, outb.first);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1- 4.*musq);
    double me2 = b*b*(5-2*b*b)*scale/120.*UnitRemoval::InvE2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(_perturbativeVertex->norm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double TFFDecayer::threeBodyME(const int , const Particle & inpart,
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
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(_tensors3, _rho3, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(_tensors3, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
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

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2,     PDT::Spin1Half,
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


  if (! (_abstractOutgoingVertex1 && _abstractOutgoingVertex2))
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr abstractOutgoingVertexF = _abstractOutgoingVertex1;
  AbstractFFVVertexPtr abstractOutgoingVertexA = _abstractOutgoingVertex2;

  if(_abstractOutgoingVertex1!=_abstractOutgoingVertex2 &&
     _abstractOutgoingVertex1->isIncoming(getParticleData(decay[ianti]->id())))
    swap (abstractOutgoingVertexF, abstractOutgoingVertexA);  
  
  if(! (inpart.dataPtr()->iColour()==PDT::Colour0)){
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;
  }

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int ifm = 0; ifm < 2; ++ifm) {
      for(unsigned int ia = 0; ia < 2; ++ia) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

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
		<< interS        .particle()->PDGName() << " in TFFDecayer::threeBodyME"
		<< Exception::runerror;

	    double gs    =  abstractOutgoingVertexF->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,_wave3[ia], interS,_tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, ifm, ia, ig) += 
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
		<< interS    .particle()->PDGName() << " in TFFDecayer::threeBodyME"
		<< Exception::runerror;
	    
	    double gs    =  abstractOutgoingVertexA->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,interS,_wavebar3[ifm],_tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, ifm, ia, ig) += 
		colourFlow[2][ix].second*diag;
	    }
	  }

	  // radiation from 4 point vertex
	  if (_abstractFourPointVertex){
	    double gs    = _abstractFourPointVertex->strongCoupling(scale);
	    Complex diag = _abstractFourPointVertex->evaluate(scale, _wave3[ia], _wavebar3[ifm],
							      _gluon[2*ig], _tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, ifm, ia, ig) += 
		colourFlow[3][ix].second*diag;
	    }
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

