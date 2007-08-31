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
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

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
    iferm = 0;
    ianti = 1;
  }

  vector<SpinorWaveFunction> awave;
  vector<SpinorBarWaveFunction> fwave;
  SpinorWaveFunction(awave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave,decay[iferm],outgoing,true,vertex);
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti){
	newme(0, ia, ifm) = _theFFSPtr->evaluate(scale,awave[ia],
					     fwave[ifm],inwave);
      }
      else {
	newme(0, ifm, ia) = _theFFSPtr->evaluate(scale,awave[ia],
					     fwave[ifm],inwave);

      }
    }
  }
  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if(decay[0]->coloured() && decay[1]->coloured())
    output *= 3.;  
  if( decay[0]->id() == decay[1]->id() )
    output *= 0.5;
  colourConnections(inpart, decay);
  return output;
}

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  _theFFSPtr->setCoupling(sqr(inpart.second), outb.first, outa.first,
			  inpart.first, 3);
  double mu1(outa.second/inpart.second),
    mu2(outb.second/inpart.second);
  double c2 = norm(_theFFSPtr->getNorm());
  Complex al(_theFFSPtr->getLeft()), ar(_theFFSPtr->getRight());
  double me2 = -c2*( (norm(al) + norm(ar))*( sqr(mu1) + sqr(mu2) - 1.)
		     + 2.*(ar*conj(al) + al*conj(ar)).real()*mu1*mu2 );
  Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
				      outb.second);
  Energy output = me2*pcm/(8*Constants::pi);
  int cola(outa.first->iColour()), colb(outa.first->iColour());
  if( abs(cola) == 3 && abs(colb) == 3)
    output *= 3.;
  if( outa.first->id() == outb.first->id() )
    output *= 0.5;
  return output;
}

