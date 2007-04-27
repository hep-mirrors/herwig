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
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
 
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
  double output=(newme.contract(rhoin)).real()/scale;
  if(decay[0]->coloured()){
    output*=3.;
  }
  colourConnections(inpart, decay);
  return output;
}

double VFFDecayer::partialWidth(const PDPtr inpart, const PDPtr anti,
				const PDPtr ferm) const {
  //Use analytic expression for partial width
  double mu1 = ferm->mass()/inpart->mass();
  double mu2 = anti->mass()/inpart->mass();
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),anti->mass(),
				      ferm->mass());
  Energy2 q2(inpart->mass()*inpart->mass());
  //coupling may not have been calculated yet
  _theFFVPtr->setCoupling(q2,anti,ferm,inpart);
  Complex norm = (_theFFVPtr->getNorm()*_theFFVPtr->getNorm());
  double cl(_theFFVPtr->getLeft().real()),cr(_theFFVPtr->getRight().real());
  double matrixElement2 = (cl*cl+cr*cr)*((mu1*mu1 + mu2*mu2 )*(mu1*mu1 + mu2*mu2) + mu2*mu2 + mu1*mu1 + 2);
  matrixElement2 += -12.*cl*cr*mu1*mu2;
  matrixElement2 *= norm.real()/3.;
  double output = matrixElement2*pcm/(8*Constants::pi);
  if(ferm->coloured())
    output *= 3.;
  return output;
}


