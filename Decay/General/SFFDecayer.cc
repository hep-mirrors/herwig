// -*- C++ -*-
//
// SFFDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
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

IBPtr SFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SFFDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<FFSVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractFFSVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<SFFDecayer> SFFDecayer::initSFFDecayer;
// Definition of the static class description member.

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

double SFFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const {
  // work out which is the fermion and antifermion
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    ME(DecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half));
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
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
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti){
	ME()(0, ia, ifm) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);
      }
      else {
	ME()(0, ifm, ia) = _abstractVertex->evaluate(scale,_wave[ia],
						     _wavebar[ifm],_swave);
	
      }
    }
  }
  double output = (ME().contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    _perturbativeVertex->setCoupling(sqr(inpart.second), outb.first, outa.first,
				     in);
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    double c2 = norm(_perturbativeVertex->norm());
    Complex al(_perturbativeVertex->left()), ar(_perturbativeVertex->right());
    double me2 = -c2*( (norm(al) + norm(ar))*( sqr(mu1) + sqr(mu2) - 1.)
		       + 2.*(ar*conj(al) + al*conj(ar)).real()*mu1*mu2 );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/(8*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

