// -*- C++ -*-
//
// PartonSplitter.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include "ClusterHadronizationHandler.h"

using namespace Herwig;

IBPtr PartonSplitter::clone() const {
  return new_ptr(*this);
}

IBPtr PartonSplitter::fullclone() const {
  return new_ptr(*this);
}

void PartonSplitter::persistentOutput(PersistentOStream & os) const {
  os << _quarkSelector << ounit(_gluonDistance,femtometer);
}

void PartonSplitter::persistentInput(PersistentIStream & is, int) {
  is >> _quarkSelector >> iunit(_gluonDistance,femtometer);
}

DescribeClass<PartonSplitter,Interfaced> 
describePartonSplitter("Herwig::PartonSplitter","");

void PartonSplitter::Init() {

  static ClassDocumentation<PartonSplitter> documentation
    ("This class is reponsible of the nonperturbative splitting of partons");
 
}

void PartonSplitter::split(PVector & tagged) {
  // set the gluon c tau once and for all
  static bool first = true;
  if(first) {
    _gluonDistance = hbarc*getParticleData(ParticleID::g)->constituentMass()/
      ClusterHadronizationHandler::currentHandler()->minVirtuality2();
    first = false;
  }
  PVector newtag;
  Energy2 Q02 = 0.99*sqr(getParticleData(ParticleID::g)->constituentMass());
  // Loop over all of the particles in the event.
  for(PVector::iterator pit = tagged.begin(); pit!=tagged.end(); ++pit) {
    // only considering gluons so add other particles to list of particles
    if( (**pit).data().id() != ParticleID::g ) {
      newtag.push_back(*pit);
      continue;
    }
    // should not have been called for massless or space-like gluons
    if((**pit).momentum().m2() <= 0.0*sqr(MeV) ) {
      throw Exception()
	<< "Spacelike or massless gluon m2= " << (**pit).momentum().m2()/GeV2
	<< "GeV2 in PartonSplitter::split()"
	<< Exception::eventerror;
    }
    // time like gluon gets split
    PPtr ptrQ = PPtr();
    PPtr ptrQbar = PPtr();
    splitTimeLikeGluon(*pit,ptrQ,ptrQbar);
    ptrQ->scale(Q02);
    ptrQbar->scale(Q02);

    (*pit)->colourLine()->addColoured(ptrQ);
    (*pit)->addChild(ptrQ);
    newtag.push_back(ptrQ);

    (*pit)->antiColourLine()->addAntiColoured(ptrQbar);
    (*pit)->addChild(ptrQbar);
    newtag.push_back(ptrQbar);
    
    // set the life length of gluon
    Length distance = UseRandom::rndExp(_gluonDistance);
    (**pit).setLifeLength((distance/(**pit).mass())*(**pit).momentum());
    // assume quarks same position as gluon
    ptrQ   ->setVertex((**pit).decayVertex());
    ptrQ   ->setLifeLength(Lorentz5Distance());
    ptrQbar->setVertex((**pit).decayVertex());
    ptrQbar->setLifeLength(Lorentz5Distance());
  }
  swap(tagged,newtag);
}

void PartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon, 
					PPtr & ptrQ, 
					PPtr & ptrQbar){
  // select the quark flavour
  tPDPtr quark = _quarkSelector.select(UseRandom::rnd());
  // Solve the kinematics of the two body decay  G --> Q + Qbar
  Lorentz5Momentum momentumQ;
  Lorentz5Momentum momentumQbar;
  double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
  using Constants::pi;
  double phiStar = UseRandom::rnd( -pi , pi );
  Energy constituentQmass = quark->constituentMass();

  if (ptrGluon->momentum().m() < 2.0*constituentQmass) {
    throw Exception() << "Impossible Kinematics in PartonSplitter::splitTimeLikeGluon()" 
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(ptrGluon->momentum(), constituentQmass, 
			   constituentQmass, cosThetaStar, phiStar, momentumQ, 
			   momentumQbar );
  // Create quark and anti-quark particles of the chosen flavour 
  // and set they 5-momentum (the mass is the constituent one).
  ptrQ    = new_ptr(Particle(quark      ));
  ptrQbar = new_ptr(Particle(quark->CC()));
  ptrQ    ->set5Momentum( momentumQ    );
  ptrQbar ->set5Momentum( momentumQbar );
}

void PartonSplitter::doinit() {
  Interfaced::doinit();
  // calculate the probabilties for the gluon to branch into each quark type
  // based on the available phase-space, as in fortran.
  Energy mg=getParticleData(ParticleID::g)->constituentMass();
  for( int ix=1; ix<6; ++ix ) {
    PDPtr quark = getParticleData(ix);
    Energy pcm = Kinematics::pstarTwoBodyDecay(mg,quark->constituentMass(),
					       quark->constituentMass());
    if(pcm>ZERO) _quarkSelector.insert(pcm/GeV,quark);
  }
  if(_quarkSelector.empty()) 
    throw InitException() << "At least one quark must have constituent mass less "
			  << "then the constituent mass of the gluon in "
			  << "PartonSplitter::doinit()" << Exception::runerror;
}
