// -*- C++ -*-
//
// TSSDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TSSDecayer class.
//

#include "TSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void TSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<TSSDecayer> TSSDecayer::initTSSDecayer;
// Definition of the static class description member.

void TSSDecayer::Init() {

  static ClassDocumentation<TSSDecayer> documentation
    ("This class implements the decay of a tensor particle into "
     "2 scalars.");

}

double TSSDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin2);
  rhoin.average();
  vector<LorentzTensor<double> > in;
  TensorWaveFunction(in,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  ScalarWaveFunction sca1(decay[0],outgoing,true,vertex);
  ScalarWaveFunction sca2(decay[1],outgoing,true,vertex);
  Energy2 scale(sqr(inpart.mass()));
  DecayMatrixElement newme(PDT::Spin2,PDT::Spin0,PDT::Spin0);
  for(unsigned int thel=0;thel<5;++thel) {
    TensorWaveFunction inwave(inpart.momentum(),
			      inpart.dataPtr(),
			      in[thel].xx(),in[thel].xy(),in[thel].xz(),
			      in[thel].xt(),in[thel].yx(),in[thel].yy(),
			      in[thel].yz(),in[thel].yt(),in[thel].zx(),
			      in[thel].zy(),in[thel].zz(),in[thel].zt(),
			      in[thel].tx(),in[thel].ty(),in[thel].tz(),
			      in[thel].tt());
    newme(thel,0,0) =_abstractVertex->evaluate(scale,sca1,sca2,inwave); 
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}


Energy TSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    _perturbativeVertex->setCoupling(scale, outa.first, outb.first, inpart.first);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1. - 4.*musq);
    double me2 = scale*pow(b,4)/120*UnitRemoval::InvE2;
    Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
					outb.second);
    Energy output = norm(_perturbativeVertex->getNorm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}
