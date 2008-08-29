// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

FFSDecayer::FFSDecayer() {
  addToSearchList(0);
  addToSearchList(1);
}

IBPtr FFSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFSDecayer::doinit() throw(InitException) {
  _perturbativeVertex = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _perturbativeVertex << _abstractVertex;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _perturbativeVertex >> _abstractVertex;
}

ClassDescription<FFSDecayer> FFSDecayer::initFFSDecayer;
// Definition of the static class description member.

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

double FFSDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1Half);
  
  vector<SpinorWaveFunction> wave;
  vector<SpinorBarWaveFunction> barWave;
  unsigned int iferm(0),iscal(1);
  if(decay[0]->data().iSpin() != PDT::Spin1Half) swap(iferm,iscal);
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
    if(wave[0].wave().Type() != u_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	wave[ix].conjugate();
      }
    }
  }
  else {
    SpinorBarWaveFunction(barWave,rhoin,const_ptr_cast<tPPtr>(&inpart),
			  incoming,true,vertex);
    SpinorWaveFunction(wave,decay[iferm],outgoing,true,vertex);
    if(barWave[0].wave().Type() != v_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	barWave[ix].conjugate();
      }
    }
  }
  ScalarWaveFunction scal(decay[iscal],outgoing,true,vertex);
  DecayMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) newme(if1, if2, 0) = 
	_perturbativeVertex->evaluate(scale,wave[if1],barWave[if2],scal);
      else     newme(if2, if1, 0) = 
	_abstractVertex    ->evaluate(scale,wave[if1],barWave[if2],scal);
    }
  }
  
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    double mu1(0.),mu2(0.);
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first,
				       outa.first, outb.first,1);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first,
				       outb.first, outa.first,1);
      
    }
    double c2 = norm(_perturbativeVertex->getNorm());
    Complex cl = _perturbativeVertex->getLeft();
    Complex cr = _perturbativeVertex->getRight();
    double me2 = c2*( (norm(cl) + norm(cr))*(1. + sqr(mu1) - sqr(mu2))
		      + 2.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/16./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
