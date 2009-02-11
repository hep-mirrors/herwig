// -*- C++ -*-
//
// VFFDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void VFFDecayer::doinit() throw(InitException) {
  _perturbativeVertex = dynamic_ptr_cast<FFVVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFVVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void VFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void VFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<VFFDecayer> VFFDecayer::initVFFDecayer;
// Definition of the static class description member.

void VFFDecayer::Init() {

  static ClassDocumentation<VFFDecayer> documentation
    ("The VFFDecayer implements the matrix element for the"
     " decay of a vector to fermion-antifermion pair");

}

double VFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    ME(DecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half));
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  // compute the matrix element
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int ifm = 0; ifm < 2; ++ifm) { //loop over fermion helicities
    for(unsigned int ia = 0; ia < 2; ++ia) {// loop over antifermion helicities
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {//loop over vector helicities
	if(iferm > ianti) {
	  ME()(vhel, ia, ifm) = 
	    _abstractVertex->evaluate(scale,_wave[ia],
				      _wavebar[ifm],_vectors[vhel]);
	}
	else
	  ME()(vhel,ifm,ia)=
	    _abstractVertex->evaluate(scale,_wave[ia],
				      _wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    double mu1(outa.second/inpart.second), mu2(outb.second/inpart.second);
    _perturbativeVertex->setCoupling(sqr(inpart.second), outa.first, outb.first,
				     inpart.first);
    Complex cl(_perturbativeVertex->getLeft()), cr(_perturbativeVertex->getRight());
    double me2 = (norm(cl) + norm(cr))*( sqr(sqr(mu1) - sqr(mu2)) 
					 + sqr(mu1) + sqr(mu2) - 2.)
      - 6.*(cl*conj(cr) + cr*conj(cl)).real()*mu1*mu2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(_perturbativeVertex->getNorm())*me2*pcm / 
      (24.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
