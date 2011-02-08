// -*- C++ -*-
//
// SSSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSDecayer class.
//

#include "SSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SSSDecayer::SSSDecayer() {
  addToSearchList(0);
  addToSearchList(1);
  addToSearchList(2);
}

IBPtr SSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void SSSDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<SSSVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractSSSVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void SSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void SSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<SSSDecayer> SSSDecayer::initSSSDecayer;
// Definition of the static class description member.

void SSSDecayer::Init() {

  static ClassDocumentation<SSSDecayer> documentation
    ("This class implements the decay of a scalar to 2 scalars.");

}

double SSSDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    ME(DecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0));
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      ScalarWaveFunction::
	constructSpinInfo(decay[ix],outgoing,true);
  }
  ScalarWaveFunction s1(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  ScalarWaveFunction s2(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  ME()(0,0,0) = _abstractVertex->evaluate(scale,s1,s2,_swave);
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(scale, in, outa.first, outb.first);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					       outb.second);
    double c2 = norm(_perturbativeVertex->norm());
    Energy pWidth = c2*pcm/8./Constants::pi/scale*UnitRemoval::E2;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
