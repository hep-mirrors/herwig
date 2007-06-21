// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TVVDecayer class.
//

#include "TVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzTensor;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

TVVDecayer::~TVVDecayer() {}

void TVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVVTPtr;
}

void TVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVVTPtr;
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
  rhoin.average();
  vector<LorentzTensor> in;
  bool massa(decay[0]->mass()==0),massb(decay[1]->mass()==0);
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
	newme(thel,v1hel,v2hel) = _theVVTPtr->evaluate(scale,vec1[v1hel],
						       vec2[v2hel],
						       inwave);
	if(massb) {++v2hel;}
      }
      if(massa) {++v1hel;}
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale;
  if(decay[0]->id() == decay[1]->id()) {
    output /= 2;
  }
  if(decay[0]->data().iColour()==PDT::Colour8 &&
     decay[1]->data().iColour()==PDT::Colour8) {
    output *= 8.;
  }
  colourConnections(inpart, decay);
  return output;
}
  
Energy TVVDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass()); 
  _theVVTPtr->setCoupling(scale,outa,outb,inpart);
  double mu2 = sqr(outa->mass()/inpart->mass());
  double b = sqrt(1 - 4.*mu2);
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  double me2;
  if(outa->mass() > 0. && outb->mass() > 0.) {
    me2 = scale*(30 - 20.*b*b + 3.*pow(b,4))/120.; 
  }
  else {
    me2 = scale/10.;
  }
  Complex norm(_theVVTPtr->getNorm()*_theVVTPtr->getNorm());
  me2 *= norm.real();
  Energy output = me2*pcm/(8.*Constants::pi);
  if(outa->id()==outb->id()) {
    output /=2;
  }
  if(outa->iColour()==PDT::Colour8 &&
     outb->iColour()==PDT::Colour8) {
    output *=8.;
  }
  return output;
}
