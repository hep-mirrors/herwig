// -*- C++ -*-
//
// IntrinsicPtGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IntrinsicPtGenerator class.
//

#include "IntrinsicPtGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Config/Constants.h"

#include "DipolePartonSplitter.h"

#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/PDF/HwRemDecayer.h"

using namespace Herwig;

IntrinsicPtGenerator::IntrinsicPtGenerator() 
  : HandlerBase(),
    theValenceIntrinsicPtScale(1.0*GeV),
    theSeaIntrinsicPtScale(1.0*GeV) {}

IntrinsicPtGenerator::~IntrinsicPtGenerator() {}

IBPtr IntrinsicPtGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr IntrinsicPtGenerator::fullclone() const {
  return new_ptr(*this);
}

SpinOneLorentzRotation IntrinsicPtGenerator::kick(PPair& in,
						  PList& intermediates) {

  if ( theValenceIntrinsicPtScale == 0.0*GeV &&
       theSeaIntrinsicPtScale == 0.0*GeV )
    return SpinOneLorentzRotation();

  assert(ShowerHandler::currentHandler());

  tHwRemDecPtr remDec = ShowerHandler::currentHandler()->remnantDecayer();

  assert(remDec);

  Lorentz5Momentum Q = in.first->momentum() + in.second->momentum();

  // first parton

  if (in.first->coloured()) {

    Axis perp;
    Energy pt = 0.*GeV;

    double phi = 2.*Constants::pi*UseRandom::rnd();

    Axis dir = in.first->momentum().vect().unit();
    perp = dir.orthogonal();
    perp.rotate(phi,dir);

    double r = sqrt(-log(1.-UseRandom::rnd()));

    if ( remDec->content().first.isValenceQuark(in.first) )
      pt = sqrt(2.) * theValenceIntrinsicPtScale * r;
    else {
      assert(in.first->id() == ParticleID::g || 
	     remDec->content().first.isSeaQuark(in.first));
      pt = sqrt(2.) * theSeaIntrinsicPtScale * r;
    }

    PPtr nin = new_ptr(Particle(in.first->dataPtr())); 

    DipolePartonSplitter::change(in.first,nin,true);

    nin->set5Momentum(Lorentz5Momentum(0.*GeV,in.first->momentum().vect() +
				       pt * perp));

    intermediates.push_back(in.first);
    in.first = nin;

  }

  // second parton

  if (in.second->coloured()) {

    Axis perp;
    Energy pt = 0.*GeV;

    double phi = 2.*Constants::pi*UseRandom::rnd();

    Axis dir = in.second->momentum().vect().unit();
    perp = dir.orthogonal();
    perp.rotate(phi,dir);

    double r = sqrt(-log(1.-UseRandom::rnd()));

    if ( remDec->content().second.isValenceQuark(in.second) )
      pt = sqrt(2.) * theValenceIntrinsicPtScale * r;
    else {
      assert(in.second->id() == ParticleID::g || 
	     remDec->content().second.isSeaQuark(in.second));
      pt = sqrt(2.) * theSeaIntrinsicPtScale * r;
    }

    PPtr nin = new_ptr(Particle(in.second->dataPtr())); 
    nin->colourInfo(new_ptr(ColourBase()));

    DipolePartonSplitter::change(in.second,nin,true);

    nin->set5Momentum(Lorentz5Momentum(0.*GeV,in.second->momentum().vect() +
				       pt * perp));

    intermediates.push_back(in.second);
    in.second = nin;

  }

  // restore mometum conservation

  Lorentz5Momentum nQ = in.first->momentum() + in.second->momentum();

  double x = Q.m()/nQ.m();

  nQ *= x;

  // work around for constructor problems
  // in Lorentz5Vector constructor, mass is not guaranteed
  // to be zero, event it was set as such before,
  // since gets recalculated in the LorentzVector
  Lorentz5Momentum scaled = x * in.first->momentum(); 
  scaled.setMass(ZERO); scaled.rescaleEnergy();
  in.first->set5Momentum(scaled);
  scaled = x * in.second->momentum(); 
  scaled.setMass(ZERO); scaled.rescaleEnergy();
  in.second->set5Momentum(scaled);

  // apparently, this is more stable than boosting directly with
  // n_Q.boostVector()-Q.boostVector()

  Boost beta1 = -Q.boostVector();
  Boost beta2 = nQ.boostVector();

  SpinOneLorentzRotation transform (beta1);
  transform.boost(beta2);

  return transform;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IntrinsicPtGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(theValenceIntrinsicPtScale,GeV) << ounit(theSeaIntrinsicPtScale,GeV);
}

void IntrinsicPtGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theValenceIntrinsicPtScale,GeV) >> iunit(theSeaIntrinsicPtScale,GeV);
}

ClassDescription<IntrinsicPtGenerator> IntrinsicPtGenerator::initIntrinsicPtGenerator;
// Definition of the static class description member.

void IntrinsicPtGenerator::Init() {

  static ClassDocumentation<IntrinsicPtGenerator> documentation
    ("IntrinsicPtGenerator generates intrinsic pt for massless "
     "incoming partons in a shower independent way.");


  static Parameter<IntrinsicPtGenerator,Energy> interfaceValenceIntrinsicPtScale
    ("ValenceIntrinsicPtScale",
     "The width of the intrinsic pt Gaussian distribution for valence partons.",
     &IntrinsicPtGenerator::theValenceIntrinsicPtScale, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<IntrinsicPtGenerator,Energy> interfaceSeaIntrinsicPtScale
    ("SeaIntrinsicPtScale",
     "The width of the intrinsic pt Gaussian distribution for sea partons.",
     &IntrinsicPtGenerator::theSeaIntrinsicPtScale, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

