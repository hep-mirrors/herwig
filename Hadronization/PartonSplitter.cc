// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonSplitter class.
//

#include "PartonSplitter.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;

void PartonSplitter::persistentOutput(PersistentOStream &) const {
}


void PartonSplitter::persistentInput(PersistentIStream &, int) {
}


ClassDescription<PartonSplitter> PartonSplitter::initPartonSplitter;
// Definition of the static class description member.


void PartonSplitter::Init() {

  static ClassDocumentation<PartonSplitter> documentation
    ("This class is reponsible of the nonperturbative splitting of partons");
 
}


tPVector PartonSplitter::split(const tPVector & tagged, tStepPtr pstep) {
  tPVector newtag;
  // Loop over all of the particles in the event.
  for(tPVector::const_iterator pit = tagged.begin(); pit!=tagged.end(); ++pit) {
    // only considering gluons so add other particles to list of particles
    if( (**pit).data().id() != ParticleID::g ) {
      newtag.push_back(*pit);
      continue;
    }
    // should not have been called for massless or space-like gluons
    if((**pit).momentum().m2() <= 0.0 ) {
      throw Exception()
	<< "Spacelike or massless gluon m2= " << (**pit).momentum().m2()/GeV2
	<< "GeV2 in PartonSplitter::split()"
	<< Exception::eventerror;
    }
    // time like gluon gets split
    PPtr ptrQ = PPtr();
    PPtr ptrQbar = PPtr();
    splitTimeLikeGluon(*pit,ptrQ,ptrQbar);
    Energy Q0 = getParticleData(ParticleID::g)->constituentMass();
    ptrQ->scale(Q0*Q0);
    ptrQbar->scale(Q0*Q0);
    pstep->addDecayProduct(*pit,ptrQ);
    pstep->addDecayProduct(*pit,ptrQbar);
    newtag.push_back(ptrQ);
    newtag.push_back(ptrQbar);
    pstep->fixColourFlow();
  }
  return newtag;
}

void PartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon, PPtr & ptrQ, PPtr & ptrQbar){

  // Choose the flavour of the quark (u or d with equal 50% probability)
  long newId = 0;
  if ( UseRandom::rndbool() ) newId = ParticleID::u;
  else                        newId = ParticleID::d;
  // Solve the kinematics of the two body decay  G --> Q + Qbar
  Lorentz5Momentum momentumQ = Lorentz5Momentum();
  Lorentz5Momentum momentumQbar = Lorentz5Momentum();
  double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
  double phiStar = UseRandom::rnd( -pi , pi );
  Energy constituentQmass = getParticleData(newId)->constituentMass();

  if (ptrGluon->momentum().m() < 2.0*constituentQmass) {
    throw Exception() << "Impossible Kinematics in PartonSplitter::splitTimeLikeGluon()" 
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(ptrGluon->momentum(), constituentQmass, 
			   constituentQmass, cosThetaStar, phiStar, momentumQ, 
			   momentumQbar ); 

  // Create quark and anti-quark particles of the chosen flavour 
  // and set they 5-momentum (the mass is the constituent one).
  ptrQ    = getParticle(newId);
  ptrQbar = getParticle(-newId);
  ptrQ->set5Momentum( momentumQ );
  ptrQbar->set5Momentum( momentumQbar );
}
