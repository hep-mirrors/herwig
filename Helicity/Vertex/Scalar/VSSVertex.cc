// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSVertex class.
//

#include "VSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
AbstractNoPIOClassDescription<VSSVertex> VSSVertex::initVSSVertex;
// Definition of the static class description member.
    
void VSSVertex::Init() {
      
static ClassDocumentation<VSSVertex> documentation
  ("The VSSVertex class is hte implementation of the"
   "vector-scalar-scalar vertex for helicity amplitude calculations."
   " all such vertices should inherit from it");
 
}

// evaluate the vertex
Complex VSSVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec,
			    const ScalarWaveFunction & sca1,
			    const ScalarWaveFunction & sca2)
{
  // get the pointers to the particle data objects
  tcPDPtr Pvec  =  vec.getParticle();
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  // calculate the coupling
  setCoupling(q2,Pvec,Psca1,Psca2);
  Complex norm=getNorm();
  // calculate the vertex
  Complex vertex = -Complex(0.,1.)*norm*
    ( vec.t()*(sca1.e() -sca2.e() )-vec.x()*(sca1.px()-sca2.px()) 
      -vec.y()*(sca1.py()-sca2.py())-vec.z()*(sca1.pz()-sca2.pz()));
  return vertex;
}

// off-shell vector
VectorWaveFunction VSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const ScalarWaveFunction & sca1,
				       const ScalarWaveFunction & sca2)
{
  // get the pointers to the particle data objects
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(sca1.px()+sca2.px(),sca1.py()+sca2.py(),
					   sca1.pz()+sca2.pz(),sca1.e() +sca2.e()); 
  // calculate the coupling
  setCoupling(q2,out,Psca1,Psca2);
  Complex norm=getNorm();
  // mass and width
  Energy mass  = out->mass();
  Energy2 mass2=mass*mass;
  // calculate the prefactor
  Energy2 p2=pout.m2();
  Complex fact=getNorm()*sca1.Wave()*sca2.Wave()*propagator(iopt,p2,out);
  // compute the vector
  Complex vec[4];
  // massive outgoing vector
  if(mass!=0.)
    {
      vec[0] = fact*(sca2.px()-sca1.px());
      vec[1] = fact*(sca2.py()-sca1.py());
      vec[2] = fact*(sca2.pz()-sca1.pz());
      vec[3] = fact*(sca2.e() -sca1.e() );
    }
  // massless outgoing vector
  else
    {
      // first the dot product for the second term
      double dot = (sca1.m2()-sca2.m2())/mass2;
      // compute the vector
      vec[0] = fact*(sca2.px()-sca1.px()+dot*pout.px());
      vec[1] = fact*(sca2.py()-sca1.py()+dot*pout.py());
      vec[2] = fact*(sca2.pz()-sca1.pz()+dot*pout.pz());
      vec[3] = fact*(sca2.e() -sca1.e() +dot*pout.e() );
    }
  return VectorWaveFunction(pout,out,vec[0],vec[1],vec[2],vec[3]);
}

// return an off-shell scalar
ScalarWaveFunction VSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const VectorWaveFunction & vec,
				       const ScalarWaveFunction & sca)
{
  // momentum of the particle
  Lorentz5Momentum pout = Lorentz5Momentum(sca.px()+vec.px(),sca.py()+vec.py(),
					   sca.pz()+vec.pz(),sca.e() +vec.e() ); 
  // get pointers to the particles
  tcPDPtr Pvec = vec.getParticle();
  tcPDPtr Psca = sca.getParticle();
  // calculate the coupling
  setCoupling(q2,out,Psca,out);
  // calculate the prefactor
  Energy2 p2 = pout.m2();
  Complex fact=getNorm()*sca.Wave()*propagator(iopt,p2,out);
  // compute the wavefunction
  fact = fact*(+vec.t()*(sca.e() +pout.e() )-vec.x()*(sca.px()+pout.px())
	       -vec.y()*(sca.py()+pout.py())-vec.z()*(sca.pz()+pout.pz()));
  return ScalarWaveFunction(pout,out,fact);
}

}
}

