// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWHVertex class.
//

#include "SSWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWWHVertex::SSWWHVertex() : theh0Wfact(0.*MeV), theH0Wfact(0.*MeV), 
			     theh0Zfact(0.*MeV), theH0Zfact(0.*MeV),
			     theCoupLast(0.*MeV), theElast(0.0),
			     theq2last(0.*MeV2), theHlast(0), 
			     theGBlast(0) {
  vector<int> first, second, third;
  //ZZh0
  first.push_back(23);
  second.push_back(23);
  third.push_back(25);
  //WWh0
  first.push_back(-24);
  second.push_back(24);
  third.push_back(25);
  //ZZH0
  first.push_back(23);
  second.push_back(23);
  third.push_back(35);
  //WWH0
  first.push_back(-24);
  second.push_back(24);
  third.push_back(35);

  setList(first, second, third);
}
			     
SSWWHVertex::~SSWWHVertex() {}

void SSWWHVertex::doinit() throw(InitException) {
  VVSVertex::doinit();
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSWWHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;

  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  Energy mz = getParticleData(ParticleID::Z0)->mass();
  double sw = sqrt(theMSSM->sin2ThetaW());
  double sinalp = sin(theMSSM->higgsMixingAngle());
  double cosalp = sqrt(1. - sqr(sinalp));
  double tanbeta = theMSSM->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt( 1. - sqr(sinbeta) );
  double sinbma = sinbeta*cosalp - cosbeta*sinalp;
  double cosbma = cosbeta*cosalp - sinbeta*sinalp;
  
  theh0Wfact = mw*sinbma/sw;
  theH0Wfact = mw*cosbma/sw;
  theh0Zfact = mz*sinbma/sw/sqrt(1. - sw*sw);
  theH0Zfact = mz*cosbma/sw/sqrt(1. - sw*sw);

  orderInGem(1);
  orderInGs(0);
}

void SSWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theh0Wfact,GeV) << ounit(theH0Wfact,GeV) 
     << ounit(theh0Zfact,GeV) << ounit(theH0Zfact,GeV);
}

void SSWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theh0Wfact,GeV) >> iunit(theH0Wfact,GeV) 
     >> iunit(theh0Zfact,GeV) >> iunit(theH0Zfact,GeV);
}

ClassDescription<SSWWHVertex> SSWWHVertex::initSSWWHVertex;
// Definition of the static class description member.

void SSWWHVertex::Init() {

  static ClassDocumentation<SSWWHVertex> documentation
    ("This is the coupling of a pair of SM gauge bosons"
     "to the higgs particles in the MSSM");

}

void SSWWHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			      tcPDPtr particle3) {
  long id1(abs(particle1->id())), id2(abs(particle2->id())), 
    id3(abs(particle3->id())), higgsID(0), bosonID(0);
  if( id1 == ParticleID::h0 || id1 == ParticleID::H0 ) {
    higgsID = id1;
    bosonID = id2;
  }
  else if( id2 == ParticleID::h0 || id2 == ParticleID::H0) {
    higgsID = id2;
    bosonID = id1;
  }
  else if( id3 == ParticleID::h0 || id3 == ParticleID::H0 ) {
    higgsID = id3;
    bosonID = id1;
  }
  else {
    throw HelicityConsistencyError() 
      << "SSWWHVertex::setCoupling - The incorrect type of higgs particle "
      << "in this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::warning;
    return;
    setNorm(0.0);
  }
  if( bosonID != ParticleID::Wplus && bosonID != ParticleID::Z0 ) {
    throw HelicityConsistencyError() 
      << "SSWWHVertex::setCoupling - The gauge boson id is incorrect " 
      << bosonID << Exception::warning;
    return;
    setNorm(0.0);
  }
  assert( higgsID == ParticleID::h0 || higgsID == ParticleID::H0 );
  assert( bosonID == ParticleID::Z0 || bosonID == ParticleID::Wplus );
  if( higgsID != theHlast || bosonID != theGBlast) {
    if( higgsID == ParticleID::h0 )
      theCoupLast = (bosonID == ParticleID::Z0) ? theh0Zfact : theh0Wfact;
    else
      theCoupLast = (bosonID == ParticleID::Z0) ? theH0Zfact : theH0Wfact;
  }
  if( q2 != theq2last ) {
    theq2last = q2;
    theElast = sqrt(4.*Constants::pi*theMSSM->alphaEM(q2));
  }

  setNorm(theElast*theCoupLast*UnitRemoval::InvE);
}
