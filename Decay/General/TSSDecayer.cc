// -*- C++ -*-
//
// TSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TSSDecayer class.
//

#include "TSSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
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


void TSSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & ,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & ,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractSSTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<SSTVertexPtr>        (vert));
  }
}


void TSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_;
}

void TSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TSSDecayer,GeneralTwoBodyDecayer>
describeHerwigTSSDecayer("Herwig::TSSDecayer", "Herwig.so");

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
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&inpart),
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
    (*ME())(thel,0,0) =0.;
    for(auto vert : vertex_)
      (*ME())(thel,0,0) += vert->evaluate(scale,sca1,sca2,tensors_[thel]); 
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}


Energy TSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first, outb.first, in);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1. - 4.*musq);
    double me2 = scale*pow(b,4)/120*UnitRemoval::InvE2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

