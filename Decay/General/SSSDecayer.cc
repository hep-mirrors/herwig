// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSDecayer class.
//

#include "SSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

void SSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theSSSPtr;
}

void SSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theSSSPtr;
}

ClassDescription<SSSDecayer> SSSDecayer::initSSSDecayer;
// Definition of the static class description member.

void SSSDecayer::Init() {

  static ClassDocumentation<SSSDecayer> documentation
    ("This class implements the decay of a scalar to 2 scalars.");

}

double SSSDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),rhoin,incoming,
			    true,vertex);
  ScalarWaveFunction s1(decay[0],outgoing,true,vertex);
  ScalarWaveFunction s2(decay[1],outgoing,true,vertex);
  Energy2 scale(inpart.mass()*inpart.mass());
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newme(0,0,0) = _theSSSPtr->evaluate(scale,s1,s2,inwave);
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if(decay[0]->id() == decay[1]->id()) {
    output /=2;
  }
  if(decay[0]->coloured() && decay[1]->coloured()) {
    output*=3.;
  }
  colourConnections(inpart, decay);
  return output;
}

Energy SSSDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  _theSSSPtr->setCoupling(scale,inpart,outa,outb);
  Complex norm = (_theSSSPtr->getNorm()*_theSSSPtr->getNorm());
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  Energy pWidth = norm.real()*pcm/8./Constants::pi/scale*UnitRemoval::E2;
  if(outa->id() == outb->id()) {
    pWidth /=2;
  }
  if(outa->coloured() && outb->coloured()) {
    pWidth *= 3.;
  }
  return pWidth;
}
