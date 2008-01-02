// -*- C++ -*-
//
// VVVDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVDecayer class.
//

#include "VVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void VVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<VVVDecayer> VVVDecayer::initVVVDecayer;
// Definition of the static class description member.

void VVVDecayer::Init() {

  static ClassDocumentation<VVVDecayer> documentation
    ("The VVVDecayer class implements the decay of a vector boson "
     "into 2 vector bosons");

}

double VVVDecayer::me2(bool vertex, const int , const Particle & inpart,
                       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1);
  rhoin.average();
  vector<VectorWaveFunction> inwave,vec1,vec2;
  VectorWaveFunction(inwave,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  VectorWaveFunction(vec1,decay[0],outgoing,true,false,vertex);
  VectorWaveFunction(vec2,decay[1],outgoing,true,false,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin1);
  unsigned int iv1,iv2,iv3;
  for(iv3=0;iv3<3;++iv3) {
    for(iv2=0;iv2<3;++iv2) {
      for(iv1=0;iv1<3;++iv1) {
	newme(iv1,iv2,iv3) = _abstractVertex->evaluate(scale,vec1[iv2],
						       vec2[iv3],inwave[iv1]);
      }
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // make the colour connections
  colourConnections(inpart, decay);
  // return the answer
  return output;
}

Energy VVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first,
				     outa.first, outb.first);
    double mu1(outa.second/inpart.second), mu1sq(sqr(mu1)),
      mu2(outb.second/inpart.second), mu2sq(sqr(mu2));
    double me2 = 
      (mu1 - mu2 - 1.)*(mu1 - mu2 + 1.)*(mu1 + mu2 - 1.)*(mu1 + mu2 + 1.)
      * (sqr(mu1sq) + sqr(mu2sq) + 10.*(mu1sq*mu2sq + mu1sq + mu2sq) + 1.)
      /4./mu1sq/mu2sq;
    Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
					outb.second);
    Energy pWidth = norm(_perturbativeVertex->getNorm())*me2*pcm/8./Constants::pi;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
