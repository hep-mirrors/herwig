// -*- C++ -*-
//
// TSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TSSDecayer class.
//

#include "TSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void TSSDecayer::doinit() {
  GeneralTwoBodyDecayer::doinit();
  _perturbativeVertex = dynamic_ptr_cast<SSTVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractSSTVertexPtr>(getVertex());
}

void TSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void TSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<TSSDecayer> TSSDecayer::initTSSDecayer;
// Definition of the static class description member.

void TSSDecayer::Init() {

  static ClassDocumentation<TSSDecayer> documentation
    ("This class implements the decay of a tensor particle into "
     "2 scalars.");

}

double TSSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
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
      ScalarWaveFunction::
	constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  ScalarWaveFunction sca1(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  ScalarWaveFunction sca2(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int thel=0;thel<5;++thel) {
    (*ME())(thel,0,0) =_abstractVertex->evaluate(scale,sca1,sca2,_tensors[thel]); 
  }
  double output = (ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}


Energy TSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(scale, outa.first, outb.first, in);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1. - 4.*musq);
    double me2 = scale*pow(b,4)/120*UnitRemoval::InvE2;
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

