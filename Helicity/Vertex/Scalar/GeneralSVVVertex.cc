// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralSVVVertex class.
//

#include "GeneralSVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

// #ifdef ThePEG_TEMPLATES_IN_CC_FILE
// // #include "GeneralSVVVertex.tcc"
// #endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/EpsFunction.h"

namespace Herwig {
  namespace Helicity {
    using namespace ThePEG;
    
    GeneralSVVVertex::~GeneralSVVVertex() {}
    
    void GeneralSVVVertex::persistentOutput(PersistentOStream & os) const {}

    void GeneralSVVVertex::persistentInput(PersistentIStream & is, int) {}
    
    AbstractClassDescription<GeneralSVVVertex> GeneralSVVVertex::initGeneralSVVVertex;
    // Definition of the static class description member.
    
    void GeneralSVVVertex::Init() {
      static ClassDocumentation<GeneralSVVVertex> documentation
	("Implementation of the SVV vertex. ");
    }
    
    Complex GeneralSVVVertex::evaluate(Energy2 q2, 
				       const ScalarWaveFunction & sca,
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2) {
      Lorentz5Momentum pSca = sca.getMomentum();
      Lorentz5Momentum pVec1 = vec1.getMomentum();
      Lorentz5Momentum pVec2 = vec2.getMomentum();
      if(kinematics()) 	{
	  calculateKinematics(pSca,pVec1,pVec2);
	}
      setCoupling(q2,sca.getParticle(),vec1.getParticle(),
		  vec2.getParticle());
      complex<Energy> e1p1(vec1.Wave()*pVec1),e1p2(vec1.Wave()*pVec2);
      complex<Energy> e2p1(vec2.Wave()*pVec1),e2p2(vec2.Wave()*pVec2);
      LorentzPolarizationVector eps;
      eps = EpsFunction::product(vec1.Wave(),vec2.Wave(),
				 pVec2);
      complex<Energy> p1Ep2 =eps*pVec1;
      Complex output=_a00*(vec1.Wave()*vec2.Wave())
	+ _a11*e1p1*e2p1 + _a12*e1p1*e2p2
	+ _a21*e1p2*e2p1 + _a22*e1p2*e2p2
	+ _aEp*p1Ep2;
      return Complex(0.,1.)*getNorm()*sca.Wave()*output;
    }
  }
}
