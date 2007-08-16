// -*- C++ -*-
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
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzTensor;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

TSSDecayer::~TSSDecayer() {}

void TSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theSSTPtr;
}

void TSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theSSTPtr;
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
  Energy2 scale(inpart.scale());
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
    newme(thel,0,0) =_theSSTPtr->evaluate(scale,sca1,sca2,inwave); 
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if(decay[0]->id() == decay[1]->id()) {
    output /=2.;
  }
  return output;
}


Energy TSSDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  _theSSTPtr->setCoupling(scale,outa,outb,inpart);
  double musq = sqr(outa->mass()/inpart->mass());
  double b = sqrt(1.-4.*musq);
  double me2 = scale*pow(b,4)/120*UnitRemoval::InvE2;
  Energy pcm = Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				      outb->mass());
  Energy output = norm(_theSSTPtr->getNorm())*me2*pcm/(8.*Constants::pi);
  if(outa->id() == outb->id()) {
    output /= 2;
  }
  return output;
}
