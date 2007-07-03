// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

SFFDecayer::~SFFDecayer() {}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFSPtr;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFSPtr;
}

ClassDescription<SFFDecayer> SFFDecayer::initSFFDecayer;
// Definition of the static class description member.

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

double SFFDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  ScalarWaveFunction inwave =
    ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),rhoin,
		       incoming,true,vertex);
  // work out which is the fermion and antifermion
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) {
    iferm=0;
    ianti=1;
  }

  vector<SpinorWaveFunction> awave;
  vector<SpinorBarWaveFunction> fwave;
  SpinorWaveFunction(awave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave,decay[iferm],outgoing,true,vertex);
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int ifm,ia;
  for(ifm=0;ifm<2;++ifm){
    for(ia=0;ia<2;++ia) {
      if(iferm>ianti){
	newme(0,ia,ifm)=_theFFSPtr->evaluate(scale,awave[ia],
					     fwave[ifm],inwave);
      }
      else {
	newme(0,ifm,ia)=_theFFSPtr->evaluate(scale,awave[ia],
					     fwave[ifm],inwave);

      }
    }
  }
  ME(newme);
  double output=(newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if(decay[0]->coloured() && decay[1]->coloured()) {
    output*=3.;
  }
  colourConnections(inpart, decay);
  return output;
}

Energy SFFDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  double mu1(outa->mass()/inpart->mass()),mu2(outb->mass()/inpart->mass());
  Energy2 q2(inpart->mass()*inpart->mass());
  _theFFSPtr->setCoupling(q2,outb,outa,inpart,3);
  Complex norm(_theFFSPtr->getNorm()*_theFFSPtr->getNorm());
  Complex cl = _theFFSPtr->getLeft();
  Complex cr = _theFFSPtr->getRight();
  Complex dl = conj(_theFFSPtr->getLeft());
  Complex dr = conj(_theFFSPtr->getRight());
  double matrixElement2 = -( (cl*dl + cr*dr)*(mu2*mu2 + mu1*mu1 - 1.)
			     -2.*mu1*mu2*(cr*dl + cl*dr) ).real();
  matrixElement2 *= norm.real();
  Energy pcm(Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				    outb->mass()));
  Energy output = matrixElement2*pcm/(8*Constants::pi);
  if(outa->coloured() && outb->coloured()){
    output *= 3.;
  }
  return output;
}

