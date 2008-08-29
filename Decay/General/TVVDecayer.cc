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

void TVVDecayer::doinit() throw(InitException) {
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

double TVVDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin2);
  
  vector<LorentzTensor<double> > in;
  bool massa(decay[0]->mass()==0*MeV),massb(decay[1]->mass()==0*MeV);
  TensorWaveFunction(in,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  vector<VectorWaveFunction> vec1,vec2;
  VectorWaveFunction(vec1,decay[0],outgoing,true,massa,vertex);
  VectorWaveFunction(vec2,decay[1],outgoing,true,massb,vertex);
  DecayMatrixElement newme(PDT::Spin2,PDT::Spin1,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int thel,v1hel,v2hel;
  for(thel=0;thel<5;++thel) {
    TensorWaveFunction inwave(inpart.momentum(),
			      inpart.dataPtr(),
			      in[thel].xx(),in[thel].xy(),in[thel].xz(),
			      in[thel].xt(),in[thel].yx(),in[thel].yy(),
			      in[thel].yz(),in[thel].yt(),in[thel].zx(),
			      in[thel].zy(),in[thel].zz(),in[thel].zt(),
			      in[thel].tx(),in[thel].ty(),in[thel].tz(),
			      in[thel].tt());
    for(v1hel=0;v1hel<3;++v1hel) {
      for(v2hel=0;v2hel<3;++v2hel) {
	newme(thel,v1hel,v2hel) = _abstractVertex->evaluate(scale,vec1[v1hel],
							    vec2[v2hel],
							    inwave);
	if(massb) ++v2hel;
      }
      if(massa) ++v1hel;
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}
  
Energy TVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    _perturbativeVertex->setCoupling(scale, outa.first, outb.first, inpart.first);
    double mu2 = sqr(outa.second/inpart.second);
    double b = sqrt(1 - 4.*mu2);
    Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
					outb.second);
    Energy2 me2;
    if(outa.second > 0.*MeV && outb.second > 0.*MeV)
      me2 = scale*(30 - 20.*b*b + 3.*pow(b,4))/120.; 
    else 
      me2 = scale/10.;
    
    Energy output = norm(_perturbativeVertex->getNorm())*me2*pcm
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
