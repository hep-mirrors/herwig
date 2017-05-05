// -*- C++ -*-
//
// TVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TVVDecayer class.
//

#include "TVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void TVVDecayer::doinit() {
  GeneralTwoBodyDecayer::doinit();
  _perturbativeVertex        = dynamic_ptr_cast<VVTVertexPtr>         (getVertex());
  _abstractVertex            = dynamic_ptr_cast<AbstractVVTVertexPtr> (getVertex());
  _abstractOutgoingVertex1   = dynamic_ptr_cast<AbstractVVVVertexPtr> (getOutgoingVertices()[0]);
  _abstractOutgoingVertex2   = dynamic_ptr_cast<AbstractVVVVertexPtr> (getOutgoingVertices()[1]);
  _abstractFourPointVertex   = dynamic_ptr_cast<AbstractVVVTVertexPtr>(getFourPointVertex());

}

void TVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex          << _perturbativeVertex
     << _abstractOutgoingVertex1 << _abstractOutgoingVertex2
     << _abstractFourPointVertex;
}

void TVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex          >> _perturbativeVertex
     >> _abstractOutgoingVertex1 >> _abstractOutgoingVertex2
     >> _abstractFourPointVertex;
}

ClassDescription<TVVDecayer> TVVDecayer::initTVVDecayer;
// Definition of the static class description member.

void TVVDecayer::Init() {

  static ClassDocumentation<TVVDecayer> documentation
    ("This class implements the decay of a tensor to 2 vector bosons");

}

double TVVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,_rho,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
  }
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(_tensors,const_ptr_cast<tPPtr>(&inpart),
			incoming,true,false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(_vectors[ix],decay[ix],outgoing,true,photon[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vectors[ix],decay[ix],outgoing,photon[ix]);
  Energy2 scale(sqr(inpart.mass()));
  unsigned int thel,v1hel,v2hel;
  for(thel=0;thel<5;++thel) {
    for(v1hel=0;v1hel<3;++v1hel) {
      for(v2hel=0;v2hel<3;++v2hel) {
	(*ME())(thel,v1hel,v2hel) = _abstractVertex->evaluate(scale,
							   _vectors[0][v1hel],
							   _vectors[1][v2hel],
							   _tensors[thel]);
	if(photon[1]) ++v2hel;
      }
      if(photon[0]) ++v1hel;
    }
  }
  double output = (ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}
  
Energy TVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(scale, outa.first, outb.first, in);
    double mu2 = sqr(outa.second/inpart.second);
    double b = sqrt(1 - 4.*mu2);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy2 me2;
    if(outa.second > ZERO && outb.second > ZERO)
      me2 = scale*(30 - 20.*b*b + 3.*pow(b,4))/120.; 
    else 
      me2 = scale/10.;
    
    Energy output = norm(_perturbativeVertex->norm())*me2*pcm
      /(8.*Constants::pi)*UnitRemoval::InvE2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double TVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay, MEOption meopt) {

  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()==ZERO; 

  // no emissions from massive vectors
  if (! (massless[0] && massless[1]))
    throw Exception()
      << "No dipoles available for massive vectors in TVVDecayer::threeBodyME"
      << Exception::runerror;

  int iglu(2);  
  if(meopt==Initialize) {
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(_tensors3, _rho3, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(_tensors3, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(_vectors3[ix],decay[ix   ],outgoing,true, massless[ix]);
    VectorWaveFunction::
        constructSpinInfo(_gluon       ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  if(nflow==2) cfactors[0][1]=cfactors[1][0];

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  // create wavefunctions
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vectors3[ix],decay[ix   ],outgoing,massless[ix]);
  VectorWaveFunction::
      calculateWaveFunctions(_gluon       ,decay[iglu ],outgoing,true);

  // // gauge test
  // _gluon.clear();
  // for(unsigned int ix=0;ix<3;++ix) {
  //   if(ix==1) _gluon.push_back(VectorWaveFunction());
  //   else {
  //     _gluon.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  // 				          decay[iglu ]->dataPtr(),10,
  // 					  outgoing));
  //   }
  // }

  
  // work out which vector each outgoing vertex corresponds to 
  if(_abstractOutgoingVertex1!=_abstractOutgoingVertex2 &&
     _abstractOutgoingVertex1->isIncoming(getParticleData(decay[1]->id())))
    swap(_abstractOutgoingVertex1, _abstractOutgoingVertex2);
  
  if (! (_abstractOutgoingVertex1 && _abstractOutgoingVertex2))
    throw Exception()
      << "Invalid vertices for QCD radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;

  if( !(inpart.dataPtr()->iColour()==PDT::Colour0))
    throw Exception()
      << "Invalid vertices for QCD radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;


  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
      for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

	  // radiation from first outgoing vector
	  if(decay[0]->dataPtr()->coloured()) {
	    assert(_abstractOutgoingVertex1);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[0]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction vectInter = 
	      _abstractOutgoingVertex1->evaluate(scale,3,off,_gluon[2*ig],
						 _vectors3[0][iv0],decay[0]->mass());
	  
	    if(_vectors3[0][iv0].particle()->PDGName()!=vectInter.particle()->PDGName())
	      throw Exception()
		<< _vectors3[0][iv0].particle()->PDGName() << " was changed to " 
		<< vectInter        .particle()->PDGName() << " in TVVDecayer::threeBodyME"
		<< Exception::runerror;

	    double gs    =  _abstractOutgoingVertex1->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,_vectors3[1][iv1], 
						     vectInter,_tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[1][ix].second*diag;
	    }
	  }

	  // radiation from second outgoing vector
	  if(decay[1]->dataPtr()->coloured()) {
	    assert(_abstractOutgoingVertex2);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[1]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction  vectInter = 
	      _abstractOutgoingVertex2->evaluate(scale,3,off,_vectors3[1][iv1],
						_gluon[2*ig],decay[1]->mass());
	    
	    if(_vectors3[1][iv1].particle()->PDGName()!=vectInter.particle()->PDGName())
	      throw Exception()
		<< _vectors3[1][iv1].particle()->PDGName() << " was changed to " 
		<< vectInter        .particle()->PDGName() << " in TVVDecayer::threeBodyME"
		<< Exception::runerror;
	    
	    double gs    =  _abstractOutgoingVertex2->strongCoupling(scale);
	    Complex diag = _abstractVertex->evaluate(scale,vectInter,_vectors3[0][iv0],
						     _tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[2][ix].second*diag;
	    }
	  }

	  // radiation from 4 point vertex
	  if (_abstractFourPointVertex){
	    double gs    = _abstractFourPointVertex->strongCoupling(scale);
	    Complex diag = _abstractFourPointVertex->evaluate(scale, _vectors3[0][iv0],
							      _vectors3[1][iv1],_gluon[2*ig], 
							      _tensors3[it])/gs;
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[3][ix].second*diag;
	    }
	  }
	}
	if(massless[1]) ++iv1;
      }
      if(massless[0]) ++iv0;
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

