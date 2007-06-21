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
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::u_spinortype;
using ThePEG::Helicity::v_spinortype;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

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
  double output = (newME.contract(rhoin)).real()/scale;
  colourConnections(inpart, decay);  
  return output;
}

Energy FFVDecayer::partialWidth(const PDPtr inpart, const PDPtr part1,
				const PDPtr part2) const {
  double mu1 = part1->mass()/inpart->mass();
  double mu2 = part2->mass()/inpart->mass();
  Energy2 q2(inpart->mass()*inpart->mass());
  _theFFVPtr->setCoupling(q2,inpart,part1,part2);
  Complex norm(_theFFVPtr->getNorm()*_theFFVPtr->getNorm());
  Complex cl(_theFFVPtr->getLeft()),cr(_theFFVPtr->getRight());
  double x = (cl*conj(cl)+cr*conj(cr)).real();
  double y = (cr*conj(cl) + cl*conj(cr)).real(); 
  double matrixElement2 = (-2*mu2*mu2 + mu1*mu1 + 1)*x - 6.*y*mu1;
  matrixElement2 += (mu1*mu1 - 1)*(mu1*mu1 - 1)*x/mu2/mu2;
  matrixElement2 *= norm.real()/2.;
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),part1->mass(),
				      part2->mass());
  Energy output =matrixElement2*pcm/(8.*Constants::pi);
  return output;
  
}

ClassDescription<FFVDecayer> FFVDecayer::initFFVDecayer;
// Definition of the static class description member.

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("There is no documentation for the FFVDecayer class");

}

