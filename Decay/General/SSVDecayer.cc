// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSVDecayer class.
//

#include "SSVDecayer.h"
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

SSVDecayer::~SSVDecayer() {}

void SSVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVSSPtr;
}

void SSVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVSSPtr;
}

ClassDescription<SSVDecayer> SSVDecayer::initSSVDecayer;
// Definition of the static class description member.

void SSVDecayer::Init() {

  static ClassDocumentation<SSVDecayer> documentation
    ("This implements the decay of a scalar to a vector and a scalar");

}

double SSVDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {

  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),rhoin,incoming,
			    true,vertex);
  unsigned int isc,ivec;
  if(decay[0]->dataPtr()->iSpin() == PDT::Spin0) {
    isc = 0;
    ivec = 1;
  }
  else {
    isc = 1;
    ivec = 0;
  }
  ScalarWaveFunction sca(decay[isc],outgoing,true,vertex);
  vector<VectorWaveFunction> vecWave;
  VectorWaveFunction(vecWave,decay[ivec],outgoing,true,false,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  //make sure decay matrix element is in the correct order
  double output(0.);
  if(ivec == 0) {
    DecayMatrixElement newme(PDT::Spin0,PDT::Spin1,PDT::Spin0);
    for(unsigned int ix = 0; ix < 3; ++ix)
      newme(0, ix, 0) = _theVSSPtr->evaluate(scale,vecWave[ix],sca, inwave);
    
    ME(newme);
    output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  }
  else {
    DecayMatrixElement newme(PDT::Spin0,PDT::Spin0,PDT::Spin1);
    for(unsigned int ix = 0; ix < 3; ++ix)
      newme(0, 0, ix) = _theVSSPtr->evaluate(scale,vecWave[ix],sca,inwave);
    
    ME(newme);
    output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  }
  colourConnections(inpart, decay);
  return output;
}

Energy SSVDecayer:: partialWidth(PMPair inpart, PMPair outa, 
				 PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  double mu1sq(0.),mu2sq(0.);
  if(outa.first->iSpin() == PDT::Spin0) {
    mu1sq = sqr(outa.second/inpart.second);
    mu2sq = sqr(outb.second/inpart.second);
    _theVSSPtr->setCoupling(sqr(inpart.second), outb.first, outa.first,
			    inpart.first);
  }
  else {
    mu1sq = sqr(outa.second/inpart.second);
    mu2sq = sqr(outb.second/inpart.second);
    _theVSSPtr->setCoupling(sqr(inpart.second), outa.first, outb.first,
			    inpart.first);
  }
  double me2(0.);
  if(mu2sq == 0.) 
    me2 = -2.*mu1sq - 2.;
  else
    me2 = ( sqr(mu2sq - mu1sq) - 2.*(mu2sq + mu1sq) + 1. )/mu2sq;
  Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
				      outb.second);
  Energy output = pcm*me2*norm(_theVSSPtr->getNorm())/8./Constants::pi;
  return output;
}
