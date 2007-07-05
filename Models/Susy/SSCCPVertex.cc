// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCCPVertex class.
//

#include "SSCCPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Models/Susy/SusyBase.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCCPVertex::SSCCPVertex() :_couplast(0.), _q2last() {
  vector<int> first, second, third(2, 22);
  for(unsigned int ic = 0; ic < 2; ++ic) {
    int chType(1000024);
    if(ic == 1) chType = 1000037;
    first.push_back(-chType);
    second.push_back(chType);
  }
  setList(first, second, third);
}

void SSCCPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS;
}

void SSCCPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS;
  _couplast = 0.;
  _q2last = 0.*GeV2;
}

ClassDescription<SSCCPVertex> SSCCPVertex::initSSCCPVertex;
// Definition of the static class description member.

void SSCCPVertex::Init() {

  static ClassDocumentation<SSCCPVertex> documentation
    ("This class implements the coupling of a photon to a pair of "
     "charginos.");

}

void SSCCPVertex::setCoupling(Energy2 q2, tcPDPtr part1 , tcPDPtr part2, 
			      tcPDPtr part3) {
  if(!(part1->id() == ParticleID::gamma ||  part2->id() == ParticleID::gamma ||
       part3->id() == ParticleID::gamma) )
    throw HelicityConsistencyError() << "There is no photon in the "
				     << "photon-chargino-chargino vertex! "
				     << Exception::warning;
  if(q2 != _q2last) {
    _q2last = q2;
    _couplast = -sqrt(4.*Constants::pi*_theSS->alphaEM(q2));
  }
  setNorm(_couplast);
  setLeft(1.);
  setRight(1.);
}
