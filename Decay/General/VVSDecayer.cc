// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSDecayer class.
//

#include "VVSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VVSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _abstractVertex << _perturbativeVertex;
}

void VVSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _abstractVertex >> _perturbativeVertex;
}

ClassDescription<VVSDecayer> VVSDecayer::initVVSDecayer;
// Definition of the static class description member.

void VVSDecayer::Init() {

  static ClassDocumentation<VVSDecayer> documentation
    ("The VVSDecayer class implements the decay of a vector"
     " to a vector and a scalar");

}

double VVSDecayer::me2(bool vertex, const int , const Particle & inpart,
 		       const ParticleVector & decay) const {
  bool massless = ( decay[0]->id()==ParticleID::gamma || 
		    decay[0]->id()==ParticleID::g );
  RhoDMatrix rhoin(PDT::Spin1);
  rhoin.average();
  vector<VectorWaveFunction> inwave,outwave;
  VectorWaveFunction(inwave,const_ptr_cast<tPPtr>(&inpart),incoming,true,
		     false,vertex);
  VectorWaveFunction(outwave,decay[0],outgoing,true,massless,vertex);
  ScalarWaveFunction sca(decay[1],outgoing,true,vertex);
  Energy2 scale(sqr(inpart.mass()));
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  for(unsigned int in=0;in<3;++in) {
    for(unsigned int out=0;out<3;++out) {
      if(massless&&out==1) ++out;
      newme(in,out,0) = _abstractVertex->evaluate(scale,inwave[in],outwave[out],sca);
    }
  }
  ME(newme);
  double output=(newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // make the colour connections
  colourConnections(inpart, decay);
  // return the answer
  return output;
}

Energy VVSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    _perturbativeVertex->setCoupling(sqr(inpart.second), inpart.first, outa.first,
				     outb.first);
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    double me2 = 2.+0.25*sqr(1.+mu1sq-mu2sq)/mu1sq;
    Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(_perturbativeVertex->getNorm())*me2*pcm/
      (8.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

