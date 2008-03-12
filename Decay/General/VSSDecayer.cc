// -*- C++ -*-
//
// VSSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSDecayer class.
//

#include "VSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void VSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<VSSDecayer> VSSDecayer::initVSSDecayer;
// Definition of the static class description member.

void VSSDecayer::Init() {

  static ClassDocumentation<VSSDecayer> documentation
    ("This implements the decay of a vector to 2 scalars");

}

double VSSDecayer::me2(bool vertex, const int , const Particle & inpart,
 		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1);
  rhoin.average();
  vector<VectorWaveFunction> inwave;
  VectorWaveFunction(inwave,const_ptr_cast<tPPtr>(&inpart),incoming,true,
		     false,vertex);
  ScalarWaveFunction sca1(decay[0],outgoing,true,vertex);
  ScalarWaveFunction sca2(decay[1],outgoing,true,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin0,PDT::Spin0);
  for(unsigned int ix=0;ix<3;++ix) {
    newme(ix,0,0) = _abstractVertex->evaluate(scale,inwave[ix],sca1,sca2);
  }
  ME(newme);
  double output=(newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first, outa.first,
				     outb.first);
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    double me2 = sqr(mu1sq - mu2sq) - 2.*(mu1sq + mu2sq);
    Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(_perturbativeVertex->getNorm())*me2*pcm /
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
