// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VFFDecayer class.
//

#include "VFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
 
VFFDecayer::~VFFDecayer() {}

void VFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFVPtr;
}

void VFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFVPtr;
}

ClassDescription<VFFDecayer> VFFDecayer::initVFFDecayer;
// Definition of the static class description member.

void VFFDecayer::Init() {

  static ClassDocumentation<VFFDecayer> documentation
    ("There is no documentation for the VFFDecayer class");

}

double VFFDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<VectorWaveFunction> inwave;
  VectorWaveFunction(inwave,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // construct the spinors for the outgoing particles
  unsigned int iferm(0), ianti(1);
  if(decay[0]->id() < 0) 
    swap(iferm, ianti);
  vector<SpinorWaveFunction> awave;
  vector<SpinorBarWaveFunction> fwave;
  SpinorWaveFunction(awave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave,decay[iferm],outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int ifm = 0; ifm < 2; ++ifm) { //loop over fermion helicities
    for(unsigned int ia = 0; ia < 2; ++ia) {// loop over antifermion helicities
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {//loop over vector helicities
	if(iferm > ianti)
	  newme(vhel, ia, ifm) = 
	    _theFFVPtr->evaluate(scale,awave[ia],
				fwave[ifm],inwave[vhel]);
	else
	  newme(vhel,ifm,ia)=
	    _theFFVPtr->evaluate(scale,awave[ia],
				fwave[ifm],inwave[vhel]);
      }
    }
  }
  ME(newme);
  double output=(newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if(decay[0]->coloured()){
    output*=3.;
  }
  colourConnections(inpart, decay);
  return output;
}

Energy VFFDecayer::partialWidth(const PDPtr inpart, const PDPtr part2,
				const PDPtr part3) const {
  //Use analytic expression for partial width
  double mu1(part2->mass()/inpart->mass()), 
    mu2(part3->mass()/inpart->mass());
  Energy2 q2(inpart->mass()*inpart->mass());
  _theFFVPtr->setCoupling(q2,part2, part3, inpart);
  Complex cl(_theFFVPtr->getLeft()), cr(_theFFVPtr->getRight());
  double me2 = (norm(cl) + norm(cr))*( sqr(sqr(mu1) - sqr(mu2)) 
				       + sqr(mu1) + sqr(mu2) - 2.)
    - 6.*(cl*conj(cr) + cr*conj(cl)).real()*mu1*mu2;
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),part2->mass(),
				      part3->mass());
  Energy output = -norm(_theFFVPtr->getNorm())*me2*pcm/(8*Constants::pi);
  if(part2->iColour() == PDT::Colour3 || part2->iColour() == PDT::Colour3bar)
    output *= 3.;
  if( part2->id() == part3->id() ) 
    output /= 2.;
  return output;
}


