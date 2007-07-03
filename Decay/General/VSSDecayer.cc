// -*- C++ -*-
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
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;


VSSDecayer::~VSSDecayer() {}

void VSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVSSPtr;
}

void VSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVSSPtr;
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
     newme(ix,0,0) = _theVSSPtr->evaluate(scale,inwave[ix],sca1,sca2);
   }
   ME(newme);
   double output=(newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
   if(decay[0]->id() == decay[1]->id()) {
     output /= 2.;
   }
   if((decay[0]->data()).iColour()==PDT::Colour3) {
     output *= 3.;
   }
   colourConnections(inpart, decay);
   return output;
 }

Energy VSSDecayer::partialWidth(const PDPtr inpart,const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  _theVSSPtr->setCoupling(scale,inpart,outa,outb);
  double mu1sq = sqr(outa->mass()/inpart->mass());
  double mu2sq = sqr(outb->mass()/inpart->mass());
  Complex norm2(_theVSSPtr->getNorm()*_theVSSPtr->getNorm());
  double me2 = (mu1sq*mu1sq + mu2sq*mu2sq - 2.*mu1sq*mu2sq - 2.*mu1sq
		- 2.*mu2sq + 1.);
  me2 *= norm2.real();
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  Energy output = me2*pcm/(8.*Constants::pi);
  if(outa->id() == outb->id()) {
    output /= 2.;
  }
  if(outa->iColour()==PDT::Colour3) {
    output *= 3.;
  }
  return output;
}

