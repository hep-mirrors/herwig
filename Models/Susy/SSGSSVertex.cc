// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSSVertex class.
//

#include "SSGSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSSVertex::SSGSSVertex() : _couplast(0.),_q2last() {
  vector<int> first,second,third;
  for(int ix=1000001;ix<1000007;++ix) {
    first.push_back(21);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(int ix=2000001;ix<2000007;++ix) {
    first.push_back(21);
    second.push_back(ix);
    third.push_back(-ix);
  }
  setList(first,second,third);
}

void SSGSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS;
}

void SSGSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS;
  _couplast=0.;
  _q2last=0.*GeV2;
}

ClassDescription<SSGSSVertex> SSGSSVertex::initSSGSSVertex;
// Definition of the static class description member.

void SSGSSVertex::Init() {

  static ClassDocumentation<SSGSSVertex> documentation
    ("There is no documentation for the SSGSSVertex class");

}

void SSGSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr part3) {
  unsigned int isf(0);
  if(part1->id() == ParticleID::g) {
    isf = abs(part2->id());
  }
  else if(part2->id() == ParticleID::g) {
    isf = abs(part1->id());
  }
  else {
    isf = abs(part1->id());
  }
  if((isf >= 1000001 && isf <= 1000006) || 
     (isf>=2000001 && isf <= 2000006) ) {
    if(q2 != _q2last) {
      double alphaStr = _theSS->alphaS(q2);
      _couplast = sqrt(4.*Constants::pi*alphaStr);
      _q2last = q2;
    }
    setNorm(_couplast);
}
  else {
    throw  HelicityConsistencyError() 
      << "SSGSSVertex::setCoupling() - Incorrect particle(s) in vertex. "
      << part1->id() << " " << part2->id() << " " <<  part3->id()
      << Exception::warning;
  }
}
