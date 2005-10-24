// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHGGVertex class.
//

#include "SMHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
/**
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
 #include "SMHGGVertex.tcc"
#endif
*/
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
  namespace Helicity {
    using namespace ThePEG;
		   
    SMHGGVertex::~SMHGGVertex() {}

    void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
      os << _theSM << _mw << _sw;
    }

    void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
      is >> _theSM >> _mw >> _sw;
      _couplast =0.;
      _q2last = 0.;
    }

    ClassDescription<SMHGGVertex> SMHGGVertex::initSMHGGVertex;
    // Definition of the static class description member.

    void SMHGGVertex::Init() {

      static ClassDocumentation<SMHGGVertex> documentation
	("This class implements the h->g,g vertex");

      static Switch<SMHGGVertex,int> interfaceQuarks
	("Quarks",
	 "Which quark loops to include in the coupling.",
	 &SMHGGVertex::_qopt, 0, false, false);
      static SwitchOption interfaceQuarksTop
	(interfaceQuarks,
	 "Top",
	 "Only include top loop",
	 0);
      static SwitchOption interfaceQuarksBottom
	(interfaceQuarks,
	 "Bottom",
	 "Include top and bottom",
	 1);
    }
    
    void SMHGGVertex::setCoupling(Energy q2, tcPDPtr part1,
				  tcPDPtr part2, tcPDPtr part3)
    {
      tcPDPtr t = getParticleData(ParticleID::t);
      tcPDPtr b = getParticleData(ParticleID::b);
      if(_qopt == 0) {
	type.resize(1,PDT::SpinUnknown);
	type[0] = PDT::Spin1Half;
	masses.resize(1,0.);
	masses[0] = _theSM->mass(q2,t);
	left.resize(1,0.);right.resize(1,0.);
	left[0] = _theSM->mass(q2,t)/_mw;
	right=left;
      }
      else {
	type.resize(2,PDT::SpinUnknown);
	type[0] = PDT::Spin1Half;
	type[1] = PDT::Spin1Half;
	masses.resize(2,0.);
	masses[0] = _theSM->mass(q2,t);
	masses[1] = _theSM->mass(q2,b);
	left.resize(2,0.);right.resize(2,0.);
        left[0] = _theSM->mass(q2,t)/_mw;
	left[1] = _theSM->mass(q2,b)/_mw;
	right = left; 
      }
      if(q2 != _q2last)	{
	double alphaStr = _theSM->alphaS(q2);
	double alpha = _theSM->alphaEM(q2);
	_couplast =2.*Constants::pi*alphaStr*sqrt(4*Constants::pi*alpha)/_sw;
	_q2last = q2;
      }
      setNorm(_couplast);
      //calculate tensor coefficients
      SVVLoopVertex::setCoupling(q2,part1,part2,part3);
    }
  }
}
 
