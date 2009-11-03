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

FFVDecayer::FFVDecayer() {
  addToSearchList(0);
  addToSearchList(1);
}

IBPtr FFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<FFVVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFVVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

double FFVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  // type of process
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;  
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(_wave,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wave[0].wave().Type() != u_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) _wave   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(_wavebar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(_wavebar[0].wave().Type() != v_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) _wavebar[ix].conjugate();
    }
    ME(DecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1));
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(_wave,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_wavebar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_wavebar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_wave,decay[0],outgoing,true);
    }
    VectorWaveFunction::
      constructSpinInfo(_vector,decay[1],outgoing,true,false);
  }
  Energy2 scale(sqr(inpart.mass()));
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(_wavebar,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(_wave   ,decay[0],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(_vector,decay[1],outgoing,false);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(ferm)
	  ME()(if1, if2,vhel) = 
	    _abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],_vector[vhel]);
	else
	  ME()(if2, if1, vhel) = 
	    _abstractVertex->evaluate(scale,_wave[if1],_wavebar[if2],_vector[vhel]);
      }
    }
  }
  double output=(ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
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
    Complex cl(_perturbativeVertex->left()),cr(_perturbativeVertex->right());
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
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = norm(_perturbativeVertex->norm())*me2*pcm/16./Constants::pi; 
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

