// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVDecayer class.
//

#include "FFVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::u_spinortype;
using ThePEG::Helicity::v_spinortype;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFVPtr;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFVPtr;
}

double FFVDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1Half);
  rhoin.average();
  DecayMatrixElement newME(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  vector<SpinorWaveFunction> wave;
  vector<SpinorBarWaveFunction> barWave;
  vector<VectorWaveFunction> vWave;
  unsigned int iferm,ivec;
  if(decay[0]->data().iSpin() == PDT::Spin1Half) {
    iferm = 0;
    ivec = 1;
  }
  else {
    iferm = 1;
    ivec=0;
  }
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(ferm) {
    SpinorWaveFunction(wave,rhoin,const_ptr_cast<tPPtr>(&inpart),
		       incoming,true,vertex);
    SpinorBarWaveFunction(barWave,decay[iferm],outgoing,true,vertex);
    //checking spinor types
    if(wave[0].wave().Type() != u_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix)
	wave[ix].conjugate();
    }
  }
  else {
    SpinorBarWaveFunction(barWave,rhoin,const_ptr_cast<tPPtr>(&inpart),
			  incoming,true,vertex);
    SpinorWaveFunction(wave,decay[iferm],outgoing,true,vertex);
    if(barWave[0].wave().Type() != v_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix)
	barWave[ix].conjugate();
    }
  }
  VectorWaveFunction(vWave,decay[ivec],outgoing,true,false,vertex);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(ferm)
	  newME(if1, if2,vhel) = _theFFVPtr->evaluate(scale,wave[if1],
						      barWave[if2],
						      vWave[vhel]);
	else
	  newME(if2, if1, vhel) = _theFFVPtr->evaluate(scale,wave[if1],
						      barWave[if2],
						      vWave[vhel]);
      }
    }
  }
  ME(newME);
  double output = (newME.contract(rhoin)).real()/scale*UnitRemoval::E2;
  colourConnections(inpart, decay);  
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  double mu1(0.),mu2(0.);
  if( outa.first->iSpin() == PDT::Spin1Half) {
    mu1 = outa.second/inpart.second;
    mu2 = outb.second/inpart.second;
    _theFFVPtr->setCoupling(sqr(inpart.second), inpart.first, outa.first,
			    outb.first);
  }
  else {
    mu1 = outb.second/inpart.second;
    mu2 = outa.second/inpart.second;
    _theFFVPtr->setCoupling(sqr(inpart.second),inpart.first, outb.first,
			    outa.first);
  }
  Complex cl(_theFFVPtr->getLeft()),cr(_theFFVPtr->getRight());
  double me2(0.);
  if( mu2 > 0. ) {
    me2 = (norm(cl) + norm(cr))*(1. +sqr(mu1*mu2) + 2.*sqr(mu2) 
				      - 2.*sqr(mu1) - 2.*sqr(mu2*mu2) 
				      +  sqr(mu1*mu1))
    - 6.*mu1*sqr(mu2)*(conj(cl)*cr + conj(cr)*cl).real();
    me2 /= sqr(mu2);
  }
  else
    me2 = 2.*( (norm(cl) + norm(cr))*(sqr(mu1) + 1.) 
	       - 4.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
  Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
				      outb.second);
  Energy output = norm(_theFFVPtr->getNorm())*me2*pcm/(8.*Constants::pi);
  return output;
  
}

ClassDescription<FFVDecayer> FFVDecayer::initFFVDecayer;
// Definition of the static class description member.

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("There is no documentation for the FFVDecayer class");

}

