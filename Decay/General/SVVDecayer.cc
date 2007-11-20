// -*- C++ -*-
//
// SVVDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVDecayer class.
//

#include "SVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::VVSVertexPtr;

SVVDecayer::~SVVDecayer() {}

void SVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVVSPtr;
}

void SVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVVSPtr;
}

ClassDescription<SVVDecayer> SVVDecayer::initSVVDecayer;
// Definition of the static class description member.

void SVVDecayer::Init() {

  static ClassDocumentation<SVVDecayer> documentation
    ("This implements the decay of a scalar to 2 vector bosons.");

}

double SVVDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  //vectors of wavefunctions to store all helicities
  vector<VectorWaveFunction> vOut1,vOut2;
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),rhoin,
			    incoming,true,vertex);
  VectorWaveFunction(vOut1,decay[0],outgoing,true,false,vertex);
  VectorWaveFunction(vOut2,decay[1],outgoing,true,false,vertex);
  //set-up DecayMatrix
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int iv1,iv2;
  for(iv2=0;iv2<3;++iv2) {
    for(iv1=0;iv1<3;++iv1) {
      newME(0,iv1,iv2) = _theVVSPtr->evaluate(scale,vOut1[iv1],vOut2[iv2],
					      inwave);
    }
  }
  ME(newME);
  double matrixElement2 = newME.contract(rhoin).real()/scale*UnitRemoval::E2;
  if(decay[0]->id() == decay[1]->id()){
    matrixElement2 /= 2.;
  } 
  colourConnections(inpart, decay);
  return matrixElement2;
}

Energy SVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  Energy2 scale(sqr(inpart.second));
  _theVVSPtr->setCoupling(scale, outa.first ,outb.first, 
			  inpart.first);
  double mu1sq = sqr(outa.second/inpart.second);
  double mu2sq = sqr(outb.second/inpart.second);
  double m1pm2 = mu1sq + mu2sq;
  double me2 = ( m1pm2*(m1pm2 - 2.) + 8.*mu1sq*mu2sq + 1.)/4./mu1sq/mu2sq;
  Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
				      outb.second);

  Energy output = norm(_theVVSPtr->getNorm())*me2*pcm/(8*Constants::pi)/scale
    *UnitRemoval::E2;
  if( outa.first->id() == outb.first->id() ) 
    output /= 2.;
  return output;
}
