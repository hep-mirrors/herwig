// -*- C++ -*-
//
// FFVDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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

using namespace Herwig;
using namespace ThePEG::Helicity;

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
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
  unsigned int iferm(0),ivec(1);
  if(decay[1]->data().iSpin() == PDT::Spin1Half) swap(iferm,ivec);
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
	  newME(if1, if2,vhel) = 
	    _abstractVertex->evaluate(scale,wave[if1],barWave[if2],vWave[vhel]);
	else
	  newME(if2, if1, vhel) = 
	    _abstractVertex->evaluate(scale,wave[if1],barWave[if2],vWave[vhel]);
      }
    }
  }
  ME(newME);
  double output=(newME.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),decay[1]->dataPtr());
  // make the colour connections
  colourConnections(inpart, decay);
  // return the answer
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    if( outa.first->iSpin() == PDT::Spin1Half)
      _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first,
				       outa.first, outb.first);
    else {
      swap(mu1,mu2);
      _perturbativeVertex->setCoupling(sqr(inpart.second),inpart.first,
				       outb.first,outa.first);
    }
    Complex cl(_perturbativeVertex->getLeft()),cr(_perturbativeVertex->getRight());
    double me2(0.);
    if( mu2 > 0. ) {
      me2 = (norm(cl) + norm(cr))*(1. + sqr(mu1*mu2) + sqr(mu2) 
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
    Energy output = norm(_perturbativeVertex->getNorm())*me2*pcm/16./Constants::pi; 
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer 
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

ClassDescription<FFVDecayer> FFVDecayer::initFFVDecayer;
// Definition of the static class description member.

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("There is no documentation for the FFVDecayer class");

}

