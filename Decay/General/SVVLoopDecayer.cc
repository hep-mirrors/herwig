// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVLoopDecayer class.
//

#include "SVVLoopDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

SVVLoopDecayer::~SVVLoopDecayer() {}

void SVVLoopDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theSVVPtr;
}
void SVVLoopDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theSVVPtr;
}

ClassDescription<SVVLoopDecayer> SVVLoopDecayer::initSVVLoopDecayer;
// Definition of the static class description member.

void SVVLoopDecayer::Init() {

  static ClassDocumentation<SVVLoopDecayer> documentation
    ("This class implements the loop mediated decay of a scalar to 2 "
     "vector bosons.");

}

double SVVLoopDecayer::me2(bool vertex, const int , 
			   const Particle & inpart,
			   const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  vector<VectorWaveFunction> vec1,vec2;
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),
			    rhoin,incoming,true,vertex);
  bool massa(decay[0]->mass()==0.*MeV);
  bool massb(decay[1]->mass()==0.*MeV);
  VectorWaveFunction(vec1,decay[0],outgoing,true,massa,vertex);
  VectorWaveFunction(vec2,decay[1],outgoing,true,massb,vertex);
  //Set up decay matrix
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int v1hel,v2hel;
  for(v1hel = 0;v1hel < 3;v1hel+=1) {
    if(massa && v1hel==1){++v1hel;}
    for(v2hel = 0;v2hel < 3;v2hel+=1) {
      if(massa && v2hel==1){++v2hel;}
      newME(0,v1hel,v2hel) = _theSVVPtr->evaluate(scale,inwave,
						  vec1[v1hel],vec2[v2hel]);
    }
    
  }
  ME(newME);
  double output = newME.contract(rhoin).real()/scale*UnitRemoval::E2;
  if(decay[0]->id() == decay[1]->id()) {
    output /=2;
  }
  if(decay[0]->data().iColour() == PDT::Colour8 &&
     decay[1]->data().iColour() == PDT::Colour8){
    output *= 2.;
  }
  colourConnections(inpart, decay);
  return output;
}
  
Energy SVVLoopDecayer::partialWidth(PMPair inpart, PMPair outa, 
				    PMPair outb) const {
  Lorentz5Momentum in(0.*MeV,0.*MeV,0.*MeV,inpart.second),out1,out2;
  Kinematics::twoBodyDecay(in,outa.second,outb.second,
			   Axis(0.,0.,1.),out1,out2);
  _theSVVPtr->calculateKinematics(in,out1,out2);
  Energy2 scale(sqr(inpart.second));
  _theSVVPtr->setCoupling(scale, inpart.first,outa.first, outb.first);
  //get loop coefficients
  Complex a00(_theSVVPtr->a00());Complex a11(_theSVVPtr->a11());
  Complex a12(_theSVVPtr->a12());Complex a21(_theSVVPtr->a21());
  Complex a22(_theSVVPtr->a22());Complex aEp(_theSVVPtr->aEp());
  double mu1(outa.second/inpart.second),mu1sq(mu1*mu1);
  double mu2(outb.second/inpart.second),mu2sq(mu2*mu2);
  Energy pcm = Kinematics::CMMomentum(inpart.second,outa.second,
				      outb.second);
  Complex me2 = (aEp*aEp + a11*(2.*a11-2.*a21+a22))*mu1sq*mu1sq
    +2.*a21*a21*mu2sq*mu1sq + 2.*mu1sq*(a11*(a21-a22) - aEp*aEp)
    + mu2sq*mu2sq*(aEp*aEp + a22*(a11-2.*a21+2.*a22)) + aEp*aEp
    + 2.*mu2sq*(aEp*aEp +(aEp*aEp +(a11-a21)*(a21-a22))*mu1sq +a11*a22-a21*a22)
    + a11*a22 + a12*(mu1sq + mu2sq -1)*(-2.*a11*mu1sq +(a21-2.*a22)*mu2sq 
					+a21*(mu1sq-1));
  me2 *= scale*scale*UnitRemoval::InvE4;
  me2 += 2.*a00*double(scale*UnitRemoval::InvE2)*((2.*a11 - a12 - a21)*mu1sq
					    - mu2sq*(a12+a21-2.*a22)
					    + a21 + a12) + 8.*a00*a00;
  me2 /= 2;
  Energy pWidth = norm(_theSVVPtr->getNorm())*me2.real()*pcm
    /(8.*Constants::pi)/scale*UnitRemoval::E2;
  if(outa.first->id() == outb.first->id())
    pWidth /= 2;
  if(outa.first->iColour() == PDT::Colour8 &&
     outb.first->iColour() == PDT::Colour8)
    pWidth *= 2.;
  return pWidth;
}
