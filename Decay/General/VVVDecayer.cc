// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVDecayer class.
//

#include "VVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

VVVDecayer::~VVVDecayer() {}

void VVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVVVPtr;
}

void VVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVVVPtr;
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
	newme(iv1,iv2,iv3) = _theVVVPtr->evaluate(scale,vec1[iv2],
						  vec2[iv3],inwave[iv1]);
      }
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale;
  if(decay[0]->id() == decay[1]->id()) {
    output /= 2.;
  }
  return output;
}

double VVVDecayer::partialWidth(const PDPtr inpart, const PDPtr outa,
                                const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  _theVVVPtr->setCoupling(scale,inpart,outa,outb);
  double mu1sq =sqr(outa->mass()/inpart->mass());
  double mu2sq =sqr(outb->mass()/inpart->mass());
  Complex norm = _theVVVPtr->getNorm()*_theVVVPtr->getNorm();
  double me2 = 1. + 2.*pow(mu1sq,4) + 9.*mu2sq - 9.*mu2sq*mu2sq + 
    6.*pow(mu1sq,3)*(3. + mu2sq) - 
    mu1sq*mu1sq*(27. + 64.*mu2sq + 18.*mu2sq*mu2sq) + 
    mu1sq*(7. - 32.*mu2sq - 2.*mu2sq*mu2sq + 10.*pow(mu2sq,3));
  me2 *= norm.real()/3./(4.*mu1sq*mu2sq);
  double pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  double pWidth = me2*pcm/8./Constants::pi;
  if(outa->id() == outb->id()) {
    pWidth /= 2.;
  }
  return pWidth;
}
