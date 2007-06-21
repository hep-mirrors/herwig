// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TFFDecayer class.
//

#include "TFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::LorentzTensor;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
 
TFFDecayer::~TFFDecayer() {}

void TFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFTPtr;
}

void TFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFTPtr;
}

ClassDescription<TFFDecayer> TFFDecayer::initTFFDecayer;
// Definition of the static class description member.

void TFFDecayer::Init() {

  static ClassDocumentation<TFFDecayer> documentation
    ("The TFFDecayer class implements the decay of a tensor particle "
     "to 2 fermions ");

}

double TFFDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin2);
  rhoin.average();
  vector<LorentzTensor> in;
  TensorWaveFunction(in,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  unsigned int iferm,ianti;
  if(decay[0]->id()<0) {
    ianti = 0;
    iferm = 1;
  }
  else {
    ianti = 1;
    iferm = 0;
  }
  vector<SpinorWaveFunction> wave;
  vector<SpinorBarWaveFunction> awave;
  SpinorWaveFunction(wave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(awave,decay[iferm],outgoing,true,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  DecayMatrixElement newme(PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half);
  unsigned int thel,fhel,ahel;
  for(thel=0;thel<5;++thel) {
    TensorWaveFunction inwave(inpart.momentum(),
			      inpart.dataPtr(),
			      in[thel].xx(),in[thel].xy(),in[thel].xz(),
			      in[thel].xt(),in[thel].yx(),in[thel].yy(),
			      in[thel].yz(),in[thel].yt(),in[thel].zx(),
			      in[thel].zy(),in[thel].zz(),in[thel].zt(),
			      in[thel].tx(),in[thel].ty(),in[thel].tz(),
			      in[thel].tt());
    for(fhel=0;fhel<2;++fhel) {
      for(ahel=0;ahel<2;++ahel) {
	if(iferm > ianti) {
	  newme(thel,fhel,ahel) = _theFFTPtr->evaluate(scale,wave[ahel],
						       awave[fhel],inwave);
	}
	else {
	  newme(thel,ahel,fhel) = _theFFTPtr->evaluate(scale,wave[ahel],
						       awave[fhel],inwave);
	}
      }
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale;
  if(decay[0]->coloured()) {
    output *= 3.;
  }
  colourConnections(inpart, decay);
  return output;
}

Energy TFFDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  _theFFTPtr->setCoupling(scale,inpart,outa,outb);
  double musq = sqr(outa->mass()/inpart->mass());
  double b = sqrt(1- 4.*musq);
  double me2 = b*b*(5-2*b*b)*scale/120;
  Complex norm(_theFFTPtr->getNorm()*_theFFTPtr->getNorm());
  me2 *= norm.real();
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  Energy pWidth = me2*pcm/(8.*Constants::pi);
  if(outa->coloured()) {
    pWidth *= 3.;
  } 
  return pWidth;
}
