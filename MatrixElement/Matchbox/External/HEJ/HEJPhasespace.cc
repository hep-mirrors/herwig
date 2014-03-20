// -*- C++ -*-
//
// HEJPhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HEJPhasespace class.
//

#include "HEJPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HEJPhasespace::HEJPhasespace() {}

HEJPhasespace::~HEJPhasespace() {}

IBPtr HEJPhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr HEJPhasespace::fullclone() const {
  return new_ptr(*this);
}

void HEJPhasespace::prepare(tStdXCombPtr xc, bool) {
  theLastXComb = xc;
}

double HEJPhasespace::generateKinematics(const double* r,
					 vector<Lorentz5Momentum>& momenta) {

  assert(lastXCombPtr()->hasMeta(HEJMetaKeys::Jets));

  CMultijet& jets = lastXCombPtr()->meta<CMultijet>(HEJMetaKeys::Jets);
  jets.reset();
  jets.setAPtype(mePartonData()[1]->id());
  jets.setBPtype(mePartonData()[0]->id());

  double weight = jets.GeneratePHSPpoint(momenta.size()-4, r);

  CLHEPConverter convert;

  convert(jets.firstIncomingMomentum(),momenta[1]);
  convert(jets.secondIncomingMomentum(),momenta[0]);

  size_t n = momenta.size();
  for ( size_t k = 2; k < n; ++k ) {
    convert(jets.outgoingMomenta()[k-2],momenta[n+1-k]);
  }

  double x1 = momenta[0].plus()/lastParticles().first->momentum().plus();
  double x2 = momenta[1].minus()/lastParticles().second->momentum().minus();

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  /*
  cerr << name() << " generated PS point\n"
       << "HEJ conventions:\n";
  cerr << jets.firstIncomingMomentum().x() << " "
       << jets.firstIncomingMomentum().y() << " "
       << jets.firstIncomingMomentum().z() << " "
       << jets.firstIncomingMomentum().t() << "\n"
       << jets.secondIncomingMomentum().x() << " "
       << jets.secondIncomingMomentum().y() << " "
       << jets.secondIncomingMomentum().z() << " "
       << jets.secondIncomingMomentum().t() << "\n";
  for ( size_t k = 0; k < momenta.size()-2; ++k )
    cerr << jets.outgoingMomenta()[k].x() << " "
	 << jets.outgoingMomenta()[k].y() << " "
	 << jets.outgoingMomenta()[k].z() << " "
	 << jets.outgoingMomenta()[k].t() << "\n";
  cerr << "Herwig++ conventions\n";
  for ( size_t k = 0; k < momenta.size(); ++k )
    cerr << momenta[k]/GeV << "\n";
  cerr << flush;
  */

  // ATTENTION check conventions
  return pow(lastSHat()/GeV2,-(double)(momenta.size()-4))*weight/(x1*x2);

}

int HEJPhasespace::nDim(int nFinal) const {
  return 3*nFinal - 1;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HEJPhasespace::persistentOutput(PersistentOStream &) const {
}

void HEJPhasespace::persistentInput(PersistentIStream &, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HEJPhasespace,Herwig::MatchboxPhasespace>
  describeHerwigHEJPhasespace("Herwig::HEJPhasespace", "HwHEJ.so");

void HEJPhasespace::Init() {

  static ClassDocumentation<HEJPhasespace> documentation
    ("HEJPhasespace provides an interface to HEJ's phasespace generator.");

}

