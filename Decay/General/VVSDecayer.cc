// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSDecayer class.
//

#include "VVSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VVSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVSDecayer::fullclone() const {
  return new_ptr(*this);
}

void VVSDecayer::doinit() {
  _perturbativeVertex = dynamic_ptr_cast<VVSVertexPtr>        (getVertex());
  _abstractVertex     = dynamic_ptr_cast<AbstractVVSVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

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

double VVSDecayer::me2(const int , const Particle & inpart,
 		       const ParticleVector & decay,
		       MEOption meopt) const {
  bool massless = ( decay[0]->id()==ParticleID::gamma || 
		    decay[0]->id()==ParticleID::g );
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors[0],_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::
      constructSpinInfo(_vectors[1],decay[0],outgoing,true,massless);
    ScalarWaveFunction::
      constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::
    calculateWaveFunctions(_vectors[1],decay[0],outgoing,massless);
  ScalarWaveFunction sca(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int in=0;in<3;++in) {
    for(unsigned int out=0;out<3;++out) {
      if(massless&&out==1) ++out;
      (*ME())(in,out,0) = 
	_abstractVertex->evaluate(scale,_vectors[0][in],_vectors[1][out],sca);
    }
  }
  double output=(ME()->contract(_rho)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VVSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(_perturbativeVertex) {
    Energy2 scale(sqr(inpart.second));
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outb.first->iSpin() == PDT::Spin0 )
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, 
				       outa.first, outb.first);
    else {
      _perturbativeVertex->setCoupling(sqr(inpart.second), in, 
				       outb.first, outa.first);
      swap(mu1sq, mu2sq);
    }
    double me2 = 2. + 0.25*sqr(1. + mu1sq - mu2sq)/mu1sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(_perturbativeVertex->norm())*me2*pcm/
      (24.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

