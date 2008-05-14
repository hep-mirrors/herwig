// -*- C++ -*-
//
// SSSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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

void SSSDecayer::doinit() throw(InitException) {
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

double SSSDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),rhoin,incoming,
			    true,vertex);
  ScalarWaveFunction s1(decay[0],outgoing,true,vertex);
  ScalarWaveFunction s2(decay[1],outgoing,true,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newme(0,0,0) = _abstractVertex->evaluate(scale,s1,s2,inwave);
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    _perturbativeVertex->setCoupling(scale, inpart.first, outa.first,
				     outb.first);
    Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
					outb.second);
    double c2 = norm(_perturbativeVertex->getNorm());
    Energy pWidth = c2*pcm/8./Constants::pi/scale*UnitRemoval::E2;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
