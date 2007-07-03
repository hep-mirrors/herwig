// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVDecayer class.
//

#include "SVVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::VVSVertexPtr;

SVVDecayer::~SVVDecayer() {}

void SVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVVSPtr;
}

void SVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVVSPtr;
}

ClassDescription<SVVDecayer> SVVDecayer::initSVVDecayer;
// Definition of the static class description member.

void SVVDecayer::Init() {

  static ClassDocumentation<SVVDecayer> documentation
    ("This implements the decay of a scalar to 2 vector bosons.");

}

double SVVDecayer::me2(bool vertex, const int , const Particle & inpart,
		       const ParticleVector & decay) const {
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  //vectors of wavefunctions to store all helicities
  vector<VectorWaveFunction> vOut1,vOut2;
  ScalarWaveFunction inwave(const_ptr_cast<tPPtr>(&inpart),rhoin,
			    incoming,true,vertex);
  VectorWaveFunction(vOut1,decay[0],outgoing,true,false,vertex);
  VectorWaveFunction(vOut2,decay[1],outgoing,true,false,vertex);
  //set-up DecayMatrix
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int iv1,iv2;
  for(iv2=0;iv2<3;++iv2) {
    for(iv1=0;iv1<3;++iv1) {
      newME(0,iv1,iv2) = _theVVSPtr->evaluate(scale,vOut1[iv1],vOut2[iv2],
					      inwave);
    }
  }
  ME(newME);
  double matrixElement2 = newME.contract(rhoin).real()/scale*UnitRemoval::E2;
  if(decay[0]->id() == decay[1]->id()){
    matrixElement2 /= 2.;
  } 
  colourConnections(inpart, decay);
  return matrixElement2;
}
  
Energy SVVDecayer::partialWidth(const PDPtr inpart,
				const PDPtr outa,
				const PDPtr outb) const {
  Energy2 scale(inpart->mass()*inpart->mass());
  Energy pcm(Kinematics::CMMomentum(inpart->mass(),outa->mass(),
				    outb->mass()));
  _theVVSPtr->setCoupling(scale,outa,outb,inpart);
  //get coupling
  Complex norm = _theVVSPtr->getNorm()*conj(_theVVSPtr->getNorm());
  double mu1(outa->mass()/inpart->mass()),mu1sq(mu1*mu1);
  double mu2(outb->mass()/inpart->mass()),mu2sq(mu2*mu2);
  double matrixElement2 = 2 + (mu1sq/mu2sq) - (1/mu2sq) + (0.25/mu1sq/mu2sq);
  matrixElement2 *= norm.real();
  Energy output = matrixElement2*pcm/(8*Constants::pi)/scale*UnitRemoval::E2;
  if(outa->id() == outb->id()) 
    output /= 2.;
  return output;
}

