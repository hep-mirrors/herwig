// -*- C++ -*-
//
// TVVDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzTensor.h"

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
  _perturbativeVertex = dynamic_ptr_cast<VVTVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractVVTVertexPtr>(getVertex());
}

void TVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void TVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
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
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,_rho,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
    ME(DecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1));
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
	ME()(thel,v1hel,v2hel) = _abstractVertex->evaluate(scale,
							   _vectors[0][v1hel],
							   _vectors[1][v2hel],
							   _tensors[thel]);
	if(photon[1]) ++v2hel;
      }
      if(photon[0]) ++v1hel;
    }
  }
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
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

