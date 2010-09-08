// -*- C++ -*-
//
// TFFDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void TFFDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<FFTVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFTVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void TFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void TFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
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
  if(decay[0]->id()>=0) swap(iferm,ianti);
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,_rho,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
    ME(DecayMatrixElement(PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half));
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
	  ME()(thel,fhel,ahel) = 
	    _abstractVertex->evaluate(scale,_wave[ahel],
				      _wavebar[fhel],_tensors[thel]);
	}
	else {
	  ME()(thel,ahel,fhel) = 
	    _abstractVertex->evaluate(scale,_wave[ahel],
				      _wavebar[fhel],_tensors[thel]);
	}
      }
    }
  }
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
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

