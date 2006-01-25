// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSSVertex class.
//

#include "SSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
SSSSVertex::~SSSSVertex() {}
    
void SSSSVertex::persistentOutput(PersistentOStream & os) const {}

void SSSSVertex::persistentInput(PersistentIStream & is, int) {}

AbstractClassDescription<SSSSVertex> SSSSVertex::initSSSSVertex;
// Definition of the static class description member.

void SSSSVertex::Init() {
  
  static ClassDocumentation<SSSSVertex> documentation
    ("The SSSSVertex class is the implementation"
     "of the helicity amplitude for the four scalar vertex"
     "all vertices of trhis type should inherit from it");
}

// evaluate the vertex
Complex SSSSVertex::evaluate(Energy2 q2, const ScalarWaveFunction & sca1,
			     const ScalarWaveFunction & sca2, 
			     const ScalarWaveFunction & sca3, 
			     const ScalarWaveFunction & sca4)
{
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  tcPDPtr Psca3 = sca3.getParticle();
  tcPDPtr Psca4 = sca4.getParticle();
  // calculate the coupling
  setCoupling(q2,Psca1,Psca2,Psca3,Psca4);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // return the answer
  return ii*norm*sca1.Wave()*sca2.Wave()*sca3.Wave()*sca4.Wave();
}

// off-shell scalar
ScalarWaveFunction SSSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					const ScalarWaveFunction & sca1,
					const ScalarWaveFunction & sca2,
					const ScalarWaveFunction & sca3)
{
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  tcPDPtr Psca3 = sca3.getParticle();
  // outgoing momentum 
  Lorentz5Momentum pout = Lorentz5Momentum(sca1.px()+sca2.px()+sca3.px(),
					   sca1.py()+sca2.py()+sca3.py(),
					   sca1.pz()+sca2.pz()+sca3.pz(),
					   sca1.e() +sca2.e() +sca3.e()); 
  // calculate the coupling
  setCoupling(q2,Psca1,Psca2,Psca3,out);
  Complex ii(0.,1.);
  // wavefunction
  Energy2 p2=pout.m2();
  Complex fact=-getNorm()*sca1.Wave()*sca2.Wave()*sca3.Wave()*
    propagator(iopt,p2,out);
  return ScalarWaveFunction(pout,out,fact);
}

}
}

