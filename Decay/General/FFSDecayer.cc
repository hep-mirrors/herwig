// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
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
using ThePEG::Helicity::u_spinortype;
using ThePEG::Helicity::v_spinortype;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


FFSDecayer::~FFSDecayer() {}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFSPtr;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFSPtr;
}

ClassDescription<FFSDecayer> FFSDecayer::initFFSDecayer;
// Definition of the static class description member.

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

double FFSDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin1Half);
  rhoin.average();
  vector<SpinorWaveFunction> wave;
  vector<SpinorBarWaveFunction> barWave;
  unsigned int iferm,iscal;
  if(decay[0]->data().iSpin() == PDT::Spin1Half) {
    iferm = 0;
    iscal = 1;
  }
  else {
    iferm = 1;
    iscal = 0;
  }
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(ferm) {
    SpinorWaveFunction(wave,rhoin,const_ptr_cast<tPPtr>(&inpart),
		       incoming,true,vertex);
    SpinorBarWaveFunction(barWave,decay[iferm],outgoing,true,vertex);
    if(wave[0].wave().Type() != u_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	wave[ix].conjugate();
      }
    }
  }
  else {
    SpinorBarWaveFunction(barWave,rhoin,const_ptr_cast<tPPtr>(&inpart),
			  incoming,true,vertex);
    SpinorWaveFunction(wave,decay[iferm],outgoing,true,vertex);
    if(barWave[0].wave().Type() != v_spinortype) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	barWave[ix].conjugate();
      }
    }
  }
  ScalarWaveFunction scal(decay[iscal],outgoing,true,vertex);
  DecayMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 scale(inpart.mass()*inpart.mass());
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) newme(if1, if2, 0) = _theFFSPtr->evaluate(scale,wave[if1],
						       barWave[if2],scal);
      else     newme(if2, if1, 0) = _theFFSPtr->evaluate(scale,wave[if1],
						       barWave[if2],scal);
    }
  }

  ME(newme);
  double output = (newme.contract(rhoin)).real()/scale*UnitRemoval::E2;
  if((inpart.data()).iColour()==PDT::Colour8) {
    output *= 0.5;
  }
  colourConnections(inpart, decay);
  
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return Energy();
  double mu1(0.),mu2(0.);
  if(outa.first->iSpin() == PDT::Spin1Half) {
    mu1 = outa.second/inpart.second;
    mu2 = outb.second/inpart.second;
    _theFFSPtr->setCoupling(sqr(inpart.second), inpart.first,
			    outa.first, outb.first,1);
  }
  else {
    mu1 = outb.second/inpart.second;
    mu2 = outa.second/inpart.second;
    _theFFSPtr->setCoupling(sqr(inpart.second), inpart.first,
			    outb.first, outa.first,1);

  }
  double c2 = norm(_theFFSPtr->getNorm());
  Complex cl = _theFFSPtr->getLeft();
  Complex cr = _theFFSPtr->getRight();
  double me2 = c2*( (norm(cl) + norm(cr))*(1. + sqr(mu1) - sqr(mu2))
		    + 2.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
  Energy pcm = Kinematics::CMMomentum(inpart.second, outa.second,
				       outb.second);
  Energy pWidth = me2*pcm/16./Constants::pi;
  if(inpart.first->iColour() == PDT::Colour8)
    pWidth *= 0.5;
  return pWidth;
}
