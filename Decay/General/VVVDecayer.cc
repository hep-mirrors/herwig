// -*- C++ -*-
//
// VVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVDecayer class.
//

#include "VVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void VVVDecayer::doinit() {
  perturbativeVertex_ = dynamic_ptr_cast<VVVVertexPtr>        (vertex());
  abstractVertex_     = dynamic_ptr_cast<AbstractVVVVertexPtr>(vertex());
  GeneralTwoBodyDecayer::doinit();
}

void VVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << abstractVertex_ << perturbativeVertex_;
}

void VVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> abstractVertex_ >> perturbativeVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVVDecayer,GeneralTwoBodyDecayer>
describeHerwigVVVDecayer("Herwig::VVVDecayer", "Herwig.so");

void VVVDecayer::Init() {

  static ClassDocumentation<VVVDecayer> documentation
    ("The VVVDecayer class implements the decay of a vector boson "
     "into 2 vector bosons");

}

double VVVDecayer::me2(const int , const Particle & inpart,
                       const ParticleVector & decay,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix) 
    massless[ix] = (decay[ix]->id()==ParticleID::gamma ||
		    decay[ix]->id()==ParticleID::g);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors_[ix+1],decay[ix],outgoing,true,massless[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix+1],decay[ix],outgoing,massless[ix]);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int iv3=0;iv3<3;++iv3) {
    for(unsigned int iv2=0;iv2<3;++iv2) {
      for(unsigned int iv1=0;iv1<3;++iv1) {
	(*ME())(iv1,iv2,iv3) = abstractVertex_->
	  evaluate(scale,vectors_[1][iv2],vectors_[2][iv3],vectors_[0][iv1]);
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_->setCoupling(sqr(inpart.second), in,
				     outa.first, outb.first);
    double mu1(outa.second/inpart.second), mu1sq(sqr(mu1)),
      mu2(outb.second/inpart.second), mu2sq(sqr(mu2));
    double vn = norm(perturbativeVertex_->norm());
    if(vn == ZERO || mu1sq == ZERO || mu2sq == ZERO) return ZERO;
    double me2 = 
      (mu1 - mu2 - 1.)*(mu1 - mu2 + 1.)*(mu1 + mu2 - 1.)*(mu1 + mu2 + 1.)
      * (sqr(mu1sq) + sqr(mu2sq) + 10.*(mu1sq*mu2sq + mu1sq + mu2sq) + 1.)
      /4./mu1sq/mu2sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy pWidth = vn*me2*pcm/24./Constants::pi;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
