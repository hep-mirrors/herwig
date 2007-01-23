// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSVertex class.
//

#include "VVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;
  
AbstractNoPIOClassDescription<VVSVertex> VVSVertex::initVVSVertex;
// Definition of the static class description member.
    
void VVSVertex::Init() {
  
  static ClassDocumentation<VVSVertex> documentation
    ("The VVSVertex class is the implementation of the"
     "vector-vector-scalar vertex. All such vertices should inherit"
     "from it.");
  
}

// evaluate the vertex
Complex VVSVertex::evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2, 
			    const ScalarWaveFunction & sca)
{
  // pointer to the particle data objects
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Psca  = sca.getParticle();
  // calculate the coupling
  setCoupling(q2,Pvec1,Pvec2,Psca);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // evaluate the vertex
  Complex vertex = ii*norm*sca.wave()*
    (vec1.t()*vec2.t()-vec1.x()*vec2.x()-vec1.y()*vec2.y()-vec1.z()*vec2.z());
  return vertex;
}

// evaluate an off-shell vector
VectorWaveFunction VVSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const VectorWaveFunction & vec,
				       const ScalarWaveFunction & sca)
{
  // pointer to the particle data objects
  tcPDPtr Pvec = vec.getParticle();
  tcPDPtr Psca = sca.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(vec.px()+sca.px(),vec.py()+sca.py(),
					   vec.pz()+sca.pz(),vec.e() +sca.e()); 
  // calculate the coupling
  setCoupling(q2,out,Pvec,Psca);
  // prefactor
  Energy2 p2=pout.m2();
  Energy mass = out->mass();
  Energy2 mass2=mass*mass;
  Complex fact=getNorm()*sca.wave()*propagator(iopt,p2,out);
  // evaluate the wavefunction
  Complex vect[4];
  // massless case
  if(mass==0.) {
    vect[0] = fact*vec.x();
    vect[1] = fact*vec.y();
    vect[2] = fact*vec.z();
    vect[3] = fact*vec.t();
  }
  // massive case
  else {
    Complex dot = (+vec.t()*pout.e() -vec.x()*pout.px()
		   -vec.y()*pout.py()-vec.z()*pout.pz())/mass2;
    vect[0] = fact*(vec.x()-dot*pout.px());
    vect[1] = fact*(vec.y()-dot*pout.py());
    vect[2] = fact*(vec.z()-dot*pout.pz());
    vect[3] = fact*(vec.t()-dot*pout.e());
  }
  return VectorWaveFunction(pout,out,vect[0],vect[1],vect[2],vect[3]);
}

// off-shell scalar
ScalarWaveFunction VVSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out, 
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2) {
  // pointer to the particle data objects
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(vec1.px()+vec2.px(),vec1.py()+vec2.py(),
					   vec1.pz()+vec2.pz(),vec1.e() +vec2.e()); 
  // calculate the coupling
  setCoupling(q2,Pvec1,Pvec2,out);
  // prefactor
  Energy2 p2=pout.m2();
  Complex fact=-getNorm()*propagator(iopt,p2,out);
  // evaluate the wavefunction
  Complex output = fact* 
    (vec1.t()*vec2.t()-vec1.x()*vec2.x()-vec1.y()*vec2.y()-vec1.z()*vec2.z());
  return ScalarWaveFunction(pout,out,output);
}

}
}

