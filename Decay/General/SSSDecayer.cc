// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSDecayer class.
//

#include "SSSDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

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

Energy SSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  Energy2 scale(sqr(inpart.second));
  _theSSSPtr->setCoupling(scale, inpart.first, outa.first,
			  outb.first);
  Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
				      outb.second);
  double c2 = norm(_theSSSPtr->getNorm());
  Energy pWidth = c2*pcm/8./Constants::pi/scale*UnitRemoval::E2;
  int cola(outa.first->iColour()), colb(outb.first->iColour());
  if( abs(cola) == 3 && abs(colb) == 3)
    pWidth *= 3.;
  if( outa.first->id() == outb.first->id() )
    pWidth *= 0.5;
  return pWidth;
}
