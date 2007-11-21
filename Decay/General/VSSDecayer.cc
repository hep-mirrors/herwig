// -*- C++ -*-
//
// VSSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
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
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


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

Energy VSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  _theVSSPtr->setCoupling(sqr(inpart.second), inpart.first, outa.first,
			  outb.first);
  double mu1sq = sqr(outa.second/inpart.second);
  double mu2sq = sqr(outb.second/inpart.second);
  double me2 = sqr(mu1sq - mu2sq) - 2.*(mu1sq + mu2sq);
  Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
				      outb.second);
  Energy output = -norm(_theVSSPtr->getNorm())*me2*pcm/(8.*Constants::pi);
  if(outa.first->id() == outb.first->id())
    output /= 2.;
  int cola(outa.first->iColour()), colb(outb.first->iColour());
  if( abs(cola) == 3 && abs(colb) == 3)
    output *= 3.;
  return output;
}

