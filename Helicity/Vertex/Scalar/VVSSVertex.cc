// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSSVertex class.
//

#include "VVSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
  
AbstractNoPIOClassDescription<VVSSVertex> VVSSVertex::initVVSSVertex;
// Definition of the static class description member.
    
void VVSSVertex::Init() {
      
static ClassDocumentation<VVSSVertex> documentation
  ("The VVSSVertex class is the implementation of helicity"
   "amplitude calculation of the vector-vector-scalar-scalar vertex."
   "All classes for this type of vertex should inherit from it.");
}
 
// evaluate the vertex
Complex VVSSVertex::evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			     const VectorWaveFunction & vec2, 
			     const ScalarWaveFunction & sca1, 
			     const ScalarWaveFunction & sca2)
{
  // pointer to the particle data objects
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  // calculate the coupling
  setCoupling(q2,Pvec1,Pvec2,Psca1,Psca2);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // evaluate the vertex
  Complex vertex = ii*norm*sca1.Wave()*sca2.Wave()*
    (vec1.t()*vec2.t()-vec1.x()*vec2.x()-vec1.y()*vec2.y()-vec1.z()*vec2.z());
  return vertex;
}

// evaluate an off-shell vector
VectorWaveFunction VVSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					const VectorWaveFunction & vec,
					const ScalarWaveFunction & sca1,
					const ScalarWaveFunction & sca2)
{
  // pointer to the particle data objects
  tcPDPtr Pvec  = vec.getParticle();
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(vec.px()+sca1.px()+sca2.px(),
					   vec.py()+sca1.py()+sca2.py(),
					   vec.pz()+sca1.pz()+sca2.pz(),
					   vec.e() +sca1.e() +sca2.e()); 
  // calculate the coupling
  setCoupling(q2,out,Pvec,Psca1,Psca2);
  Complex ii(0.,1.);
  // prefactor
  Energy2 p2=pout.m2();
  Energy mass = out->mass();
  Energy mass2=mass*mass;
  Complex fact=getNorm()*sca1.Wave()*sca2.Wave()*propagator(iopt,p2,out);
  // evaluate the wavefunction
  Complex vect[4];
  // massless case
  if(mass==0.)
    {
      vect[0] = fact*vec.x();
      vect[1] = fact*vec.y();
      vect[2] = fact*vec.z();
      vect[3] = fact*vec.t();
    }
  // massive case
  else
    {
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
ScalarWaveFunction VVSSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out, 
					const VectorWaveFunction & vec1,
					const VectorWaveFunction & vec2,
					const ScalarWaveFunction & sca)
{
  // pointer to the particle data objects
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Psca  = sca.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(vec1.px()+vec2.px()+sca.px(),
					   vec1.py()+vec2.py()+sca.py(),
					   vec1.pz()+vec2.pz()+sca.pz(),
					   vec1.e() +vec2.e() +sca.e()); 
  // calculate the coupling
  setCoupling(q2,Pvec1,Pvec2,out,Psca);
  // prefactor
  Energy2 p2=pout.m2();
  Complex fact=-getNorm()*sca.Wave()*propagator(iopt,p2,out);
  // evaluate the wavefunction
  Complex output = fact*
    (vec1.t()*vec2.t()-vec1.x()*vec2.x()-vec1.y()*vec2.y()-vec1.z()*vec2.z());
  return ScalarWaveFunction(pout,out,output);
}

}
}

